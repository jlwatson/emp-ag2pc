#include <emp-tool/emp-tool.h>
#include "test/single_execution.h"
#include <boost/mpi.hpp>

void parallel_online(C2PCDist *twopc, int party, bool *input, bool *output);

using namespace std;
using namespace emp;

int main(int argc, char** argv) {
    mpi::environment env;
    mpi::communicator world;

    cout << "Starting\n";
    int port, party;
    parse_party_and_port(argv, &party, &port);

    C2PCDist * twopc = nullptr;
    NetIO *io = nullptr;
    string file = "/home/romilb/ckts/test.txt";
    CircuitFile cf(file.c_str());

    // Perform the pre-processing steps in the main process
    if (world.rank() == 0) {
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

    // Define inputs
    bool in[2], out[2];
    memset(in, 0, sizeof(in));
    memset(out, 0, sizeof(out));

    if (party == ALICE)
        in[0] = true;
    else
        in[1] = true;

    // Run circuit
    cout << "Inputs ready, starting online phase.\n";
    fflush(stdout);
    auto t_parallel = clock_start();
    parallel_online(twopc, party, in, out);
    cout << "online:\t party: " << party << "\trank: " << world.rank() << "\t time: " << time_from(t_parallel) << endl;

    // Print the results in the master process
    if (world.rank() == 0) {

        string input = "";
        for (int i = 0; i < 2; ++i) {
            input += (in[i] ? "1" : "0");
        }

        string res = "";
        for (int i = 0; i < 2; ++i) {
            res += (out[i] ? "1" : "0");
        }

        cout << "party:" << party << ", input: " << input << endl;
        cout << "party:" << party << ", result: " << res << endl;

        delete io;
    }
    return 0;
}

void run_evaluation(C2PCDist_state *state, bool *input, uint8_t *mask_input, block * updated_labels){
    // Deserialize the state vectors and run the circuit:
    int gates[state->cf_num_gates*4];
    for (int i = 0; i < state->cf_num_gates*4; i++)
        gates[i] = state->gates[i];

    block labels[state->cf_num_wires];
    for (int i = 0; i < state->cf_num_wires; i++)
        labels[i] = state->labels[i];

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

    cout<<"Done deserializing.\n";
    cout << state->cf_num_gates;
    fflush(stdout);
    C2PCDist::bob_parallel_evaluate(input, mask_input, state->cf_num_gates, gates,
                                    state->party, labels,
                                    GT, GTK, GTM, GTv, state->fpreDelta);
    updated_labels = labels;
}

void parallel_online(C2PCDist *twopc, int party, bool *input, bool *output) {
    // The master sends/recvs masks with the other party, broadcasts state to the workers
    // worker runs the complete circuit, broadcasts the results and the master unmasks the output.

    mpi::environment env;
    mpi::communicator world;

    cout<<"Start parallel online phase.\n";
    fflush(stdout);
    int num_wires;
    if (world.rank() == 0) {
        num_wires = twopc->cf->num_wire;
    }
    broadcast(world, num_wires, 0);

    uint8_t *mask_input = new uint8_t[num_wires];
    memset(mask_input, 0, num_wires);

    C2PCDist_state state;

    if (world.rank() == 0) {
        // Alice sends the partial shares and bob updates labels
        twopc->send_recv_masks(input, mask_input);
        cout<<"Send recv done.\n";
        if (party == BOB) {
            // Update state on workers once send recv is done
            world.send(1, 42, mask_input, num_wires);
            state = C2PCDist_state(*twopc);
            world.send(1, 43, state);
        }
    }
    else{
        // Worker process in bob waits for state from master
        if (party == BOB) {
            if (world.rank() == 1) {
                world.recv(0, 42, mask_input, num_wires);   // 42 is a random tag, it can be any int.
                world.recv(0, 43, state);
            }
        }
    }

    // Evaluate circuit
    if(party == BOB) {
        if (world.rank() == 1) {
            // Worker process runs the evaluation
            cout<<"Starting ckt evaluation.\n";
            block * new_labels = nullptr;
            run_evaluation(&state, input, mask_input, new_labels);
            world.send(0, 44, mask_input, num_wires);
            cout<<"Ckt evaluation done.\n";
        }
        if (world.rank() == 0) {
            // Master process receives masks and unmasks to get final result
            cout<<"Master waiting for new masks.\n";
            world.recv(1, 44, mask_input, num_wires);
            cout<<"Ckt unmasking start.\n";
            twopc->bob_unmask_output(output, mask_input);
        }
    }
}
