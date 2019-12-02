#include <emp-tool/emp-tool.h>
#include "test/single_execution.h"

void parallel_online(C2PCDist &twopc, bool *input, bool *output);

using namespace std;
using namespace emp;

int main(int argc, char** argv) {

    mpi::environment env;
    mpi::communicator world;

    int port, party;
    parse_party_and_port(argv, &party, &port);

    if (world.rank() == 0) {
        NetIO *io = new NetIO(party == ALICE ? nullptr : IP, port);
        io->set_nodelay();
        string file = circuit_file_location + "/sort.txt";
        CircuitFile cf(file.c_str());
        auto t1 = clock_start();
        C2PCDist twopc(io, party, &cf);
        io->flush();
        cout << "one time:\t" << party << "\t" << time_from(t1) << endl;
        t1 = clock_start();
        twopc.function_independent();
        io->flush();
        cout << "inde:\t" << party << "\t" << time_from(t1) << endl;

        t1 = clock_start();
        twopc.function_dependent();
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

        t1 = clock_start();
        parallel_online(twopc, in, out);
        cout << "online:\t" << party << "\t" << time_from(t1) << endl;

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
    else {
        if (party==BOB){    // Work only if im bob
            cout << "Process " << world.rank() << " sitting idle.";
        }
    }

    return 0;
}

void parallel_online(C2PCDist &twopc, bool *input, bool *output) {
    twopc.online2(input, output);

}
