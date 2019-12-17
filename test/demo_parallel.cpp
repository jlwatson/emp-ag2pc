#include <emp-tool/emp-tool.h>
#include "test/single_execution.h"
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

#define MPI_MASK_INPUT_REQUEST 42
#define MPI_MASK_INPUT_RESPONSE 43
#define MPI_MASK_LABEL_RESPONSE 44

#define MPI_MASK_OUTPUT_WIRE_IDX 62
#define MPI_MASK_OUTPUT_VALUE 63

void parallel_online(C2PCDist *twopc, int party, bool *input, bool *output);

using namespace std;
using namespace emp;

int main(int argc, char** argv) {

    mpi::environment env;
    mpi::communicator world;

    int port, party;
    parse_party_and_port(argv, &party, &port);

    C2PCDist * twopc = nullptr;
    NetIO *io = nullptr;
    string prefix = "../circuits/partitioned_circuits/simple/simple";

    CircuitFile *cf = nullptr;
    if (party == ALICE || world.rank() == 0) {
        cf = new CircuitFile((prefix + ".txt").c_str()); 
    }

    // Read partitioned circuit files
    PartitionFile* partition_file = nullptr;
    MetadataFile *meta_file = nullptr;
    if (party == BOB) {
        partition_file = new PartitionFile((prefix + "-" + world.rank() + ".txt")); 
        meta_file = new MetadataFile((prefix + "-" + world.rank() + "-meta.txt"));
    }

    // Perform the pre-processing steps in the main process
    if (party == ALICE || world.rank() == 0) {

        NetIO *io = new NetIO(party == ALICE ? nullptr : IP, port);
        io->set_nodelay();
        auto t1 = clock_start();
        twopc = new C2PCDist(io, party, &cf);
        io->flush();
        cout << "one time:\t" << party << "\t" << time_from(t1) << endl;
        t1 = clock_start();
        twopc->function_independent();
        io->flush();
        cout << "inde:\t" << party << "\t" << time_from(t1) << endl;

        t1 = clock_start();
        twopc->function_dependent();
        io->flush();
        cout << "dep:\t" << party << "\t" << time_from(t1) << endl;
    }

    // Define inputs/outputs
    bool *in = calloc(cf->n1, sizeof(bool));
    bool *out = calloc(cf->n3, sizeof(bool));
    if (party == ALICE) {

        // Set input values
        // in[0] = true;

    } else { // BOB

        // Set input values
        // XXX: note that we are assuming we know which bits go with which rank process
        // if (world.rank() == n) {
        //     in[0] = true;
        // }
    }

    // Run circuit
    cout << "Inputs ready, starting online phase.\n";
    fflush(stdout);
    auto t_parallel = clock_start();
    // get all result bits back to rank 0 on both sides (Alice is just doing it locally anyway)
    parallel_online(twopc, party, in, out, partition_file, meta_file, env, world);
    cout << "online:\t" << party << "\t" << time_from(t_parallel) << endl;

    // Print the results in the master process
    if (world.rank() == 0) {
        string input = "";
        for (int i = 0; i < party == ALICE ? cf->n1 : cf->n2; ++i) {
            if (i % 32 == 0) {
                input += " ";
            }
            input += (in[i] ? "1" : "0");
        }

        string res = "";
        for (int i = 0; i < cf->n3; ++i) {
            if (i % 32 == 0) {
                res += " ";
            }
            res += (out[i] ? "1" : "0");
        }
        cout << "party: " << party << ", input:  " << input << endl;
        cout << "party: " << party << ", result: " << res << endl;

        delete io;
    }
    return 0;
}

void run_evaluation(C2PCDist *twopc, int party, bool *input, uint8_t *partition_mask_input, block *partition_labels, block *updated_labels, PartitionFile *partition_file, MetadataFile *meta_file, mpi::environment &env, mpi::communicator &world){

    // MPI get or deserialize the state vectors and run the circuit:

    // XXX: replace with partition_file->gates
    /*
    int gates[partition_file->num_partition_gate * 4];
    for (int i = 0; i < state->cf_num_gates*4; i++)
        gates[i] = state->gates[i];
    */

    // XXX: replace with partition_labels
    /*
    block labels[state->cf_num_wires];
    for (int i = 0; i < state->cf_num_wires; i++)
        labels[i] = state->labels[i];
    */

    block *GT;
    block *GTK;
    block *GTM;
    block *GTv;
    block *fpreDelta;

    if (world.rank() == 0) {
        GT = &twopc->GT;
        GTK = &twopc->GTK;
        GTM = &twopc->GTM;
        GTv = &twopc->GTv;
        fpreDelta = &twopc->fpreDelta;
    } else {
        GT = new block[partition_file->num_partition_and][4][2];
        GTK = new block[partition_file->num_partition_and][4];
        GTM = new block[partition_file->num_partition_and][4];
        GTv = new block[partition_file->num_partition_and][4];
        fpreDelta = new block;
    }

    if (world.rank() == 0) {
        for (int i = 0; i < twopc->num_ands - partition_file->num_partition_and; i++) {
            int and_gate_index;
            mpi::status s = world.recv(mpi::any_source, MPI_G_REQUEST, &and_gate_index);
            world.send(s.source(), MPI_G_GT_RESPONSE, &twopc->GT[and_gate_index], 4*2);
            world.send(s.source(), MPI_G_GTK_RESPONSE, &twopc->GTK[and_gate_index], 4);
            world.send(s.source(), MPI_G_GTM_RESPONSE, &twopc->GTM[and_gate_index], 4);
            world.send(s.source(), MPI_G_GTv_RESPONSE, &twopc->GTv[and_gate_index], 4);
        }

        for (int i = 0; i < world.size(); i++) {
            mpi::status s = world.recv(mpi::any_source, MPI_FPREDELTA_REQUEST);
            world.send(s.source(), MPI_FPREDELTA_RESPONSE, &twopc->fpreDelta);
        }

    } else { // receive responses
        for (int i = 0; i < partition_file->num_partition_gate; i++) {
            if (partition_file->gates[i * 4 + 3] == AND_GATE) {
                int global_and_gate_index = partition_file->local_and_gate_idx_to_global[i];
                world.send(0, MPI_G_REQUEST, &global_and_gate_index);
                world.recv(0, MPI_G_GT_RESPONSE, &GT[global_and_gate_index], 4*2);
                world.recv(0, MPI_G_GTK_RESPONSE, &GTK[global_and_gate_index], 4);
                world.recv(0, MPI_G_GTM_RESPONSE, &GTK[global_and_gate_index], 4);
                world.recv(0, MPI_G_GTv_RESPONSE, &GTK[global_and_gate_index], 4);
            }
        }

        world.send(0, MPI_FPREDELTA_REQUEST);
        world.recv(0, MPI_FPREDELTA_RESPONSE, &fpreDelta);
    }

    /*
    block GT[state->num_ands][4][2];
    for (int i = 0; i < state->num_ands; i++){
        for (int j = 0; j < 4; j++){
            for (int k = 0; k < 2; k++){
                GT[i][j][k] = state->GT[i*4*2 + j*2 + k];
            }
        }
    }

    block GTK[state->num_ands][4];
    for (int i = 0; i < state->num_ands; i++){
        for (int j = 0; j < 4; j++){
            GTK[i][j] = state->GTK[i*4 + j];
        }
    }

    block GTM[state->num_ands][4];
    for (int i = 0; i < state->num_ands; i++){
        for (int j = 0; j < 4; j++){
            GTM[i][j] = state->GTM[i*4 + j];
        }
    }

    bool GTv[state->num_ands][4];
    for (int i = 0; i < state->num_ands; i++){
        for (int j = 0; j < 4; j++){
            GTv[i][j] = state->GTv[i*4 + j];
        }
    }
    */

    C2PCDist::bob_parallel_evaluate(input, partition_mask_input, partition_file, meta_file, party, partition_labels, GT, GTK, GTM, GTv, fpreDelta, env, world);
    updated_labels = partition_labels;
}

void parallel_online(C2PCDist *twopc, int party, bool *input, bool *output, PartitionFile *partition_file, MetadataFile *meta_file, mpi::environment &env, mpi::communicator &world) {
    // The master sends/recvs masks with the other party, broadcasts state to the workers
    // worker runs the complete circuit, broadcasts the results and the master unmasks the output.
    cout << "Start parallel online phase." << endl;
    fflush(stdout);

    int *partition_mask_input = new uint8_t[num_partition_wire];
    memset(partition_mask_input, 0, num_partition_wire);

    int *partition_labels = new uint8_t[num_partition_wire];
    memset(partition_labels, 0, num_partition_wire);

    int *mask_input = nullptr; // only valid for rank 0

    if (world.rank() == 0) {
        int num_wires = partition_file->num_wires;
        mask_input = new uint8_t[num_wires];
        memset(mask_input, 0, num_wires);

        // Alice sends the partial shares and bob updates labels
        twopc->send_recv_masks(input, mask_input);
        if (party == BOB) {
            // Update mask input state on workers once send recv is done
            
            n_inputs_sent = 0;
            // while we haven't distributed all the input wires to other workers
            while (n_inputs_sent < (partition_file->num_alice_inputs + partition_file->num_bob_inputs) - meta_file->num_input_wires) {
                int wire_index;
                mpi::status s = world.recv(mpi::any_source, MPI_MASK_INPUT_REQUEST, &wire_index);
                world.send(s.source(), MPI_MASK_INPUT_RESPONSE, &mask_input[wire_index]);
                world.send(s.source(), MPI_MASK_LABEL_RESPONSE, &twopc->labels[wire_index]);
                n_inputs_sent++;
            }
        }
    }

    // Worker process in bob (including rank 0) waits for state from rank 0
    if (party == BOB) {
        for (auto &kvpair : partition_file->global_wire_idx_to_local) {
            int global_wire_idx = kvpair.first;
            if (global_wire_idx < partition_file->num_alice_inputs + partition_file->num_bob_inputs) {
                if (world.rank() == 0) {
                    partition_mask_input[kvpair.second] = mask_input[global_wire_idx];
                } else {
                    world.send(0, MPI_MASK_INPUT_REQUEST, &global_wire_idx);
                    world.recv(0, MPI_MASK_INPUT_RESPONSE, &partition_mask_input[kvpair.second]);
                    world.recv(0, MPI_MASK_LABEL_RESPONSE, &partition_labels[kvpair.second]);
                }
            }
        }

        // Evaluate the circuit partition
        // Worker process runs the evaluation
        block * new_labels = nullptr;
        run_evaluation(twopc, party, input, partition_mask_input, partition_labels, new_labels, partition_file, meta_file, env, world);

        // Send output bit masks to rank 0
        if (world.rank() != 0) {
            for (auto &kvpair : partition_file->global_wire_idx_to_local) {
                int global_wire_idx = kvpair.first;
                if (global_wire_idx > partition_file->num_wires - partition_file->num_outputs) {
                    world.send(0, MPI_MASK_OUTPUT_WIRE_IDX, &global_wire_idx);
                    world.send(0, MPI_MASK_OUTPUT_VALUE, &partition_mask_input[kvpair.second]);
                }
            }
        } else { // rank == 0
            for (int i = 0; i < partition_file->num_outputs - meta_file->num_output_wires; i++) {
                int wire_idx;
                mpi::source s = world.recv(mpi::any_source, MPI_MASK_OUTPUT_WIRE_IDX, &wire_idx);
                world.recv(s.source(), MPI_MASK_OUTPUT_VALUE, &mask_input[wire_idx]);
            }

            twopc->bob_unmask_output(output, mask_input);
        }
    }
}
