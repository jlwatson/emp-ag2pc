#include <emp-tool/emp-tool.h>
#include "test/single_execution.h"
#include <boost/mpi.hpp>

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
    string file = circuit_file_location + "/sort.txt";
    CircuitFile cf(file.c_str());

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

        bool in[64], out[64];
        memset(in, 0, sizeof(in));
        memset(out, 0, sizeof(out));

        if (party == ALICE) {
            /* 2 */
            in[6] = true;
        } else {
            /* 4 */
            in[34] = true;
        }

        auto t_parallel = clock_start();
        parallel_online(twopc, party, in, out);
        cout << "online:\t" << party << "\t" << time_from(t_parallel) << endl;

        string input = "";
        for (int i = 0; i < 64; ++i) {
            if (i % 32 == 0) {
                input += " ";
            }
            input += (in[i] ? "1" : "0");
        }

        string res = "";
        for (int i = 0; i < 64; ++i) {
            if (i % 32 == 0) {
                res += " ";
            }
            res += (out[i] ? "1" : "0");
        }
        cout << "party:" << party << ", input: " << input << endl;
        cout << "party:" << party << ", result: " << res << endl;

        delete io;
    }
    return 0;
}

void parallel_online(C2PCDist *twopc, int party, bool *input, bool *output) {

    mpi::environment env;
    mpi::communicator world;

    int num_wires = twopc->cf->num_wire;
    uint8_t *mask_input = new uint8_t[num_wires];
    memset(mask_input, 0, num_wires);

    if (world.rank() == 0) {
        twopc->send_recv_masks(input, mask_input);
    }

    // Evaluate circuit
    if(party == BOB) {
        if (world.rank() == 0) {
            // Read the circuit file and metadata file
            //
            twopc->bob_parallel_evaluate(input, mask_input, twopc->cf->num_gate, twopc->cf->gates,
                                         party, twopc->labels,
                                         twopc->GT, twopc->GTK, twopc->GTM, twopc->GTv, twopc->fpre->Delta);
            // Unmask
            twopc->bob_unmask_output(output, mask_input);
        }
    }
}
