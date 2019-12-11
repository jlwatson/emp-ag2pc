#include <emp-tool/emp-tool.h>
#include "test/single_execution.h"
using namespace std;
using namespace emp;

int main(int argc, char** argv) {


    int port, party;
    parse_party_and_port(argv, &party, &port);

    NetIO *io = new NetIO(party == ALICE ? nullptr : IP, port);
    io->set_nodelay();
    string file = argv[3];
    CircuitFile cf(file.c_str());
    auto t1 = clock_start();
    C2PC twopc(io, party, &cf);
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

    bool in[cf.n1], out[cf.n3];
    memset(in, 0, sizeof(in));
    memset(out, 0, sizeof(out));

    if (party == ALICE) {
        /* 2 */
        in[0] = true;
    } else {
        /* 4 */
        in[1] = true;
    }

    t1 = clock_start();
    twopc.online(in, out);
    cout << "online:\t party: " << party << "\t time: " << time_from(t1) << endl;

    string input = "";
    for (int i = 0; i < cf.n1; ++i) {
        input += (in[i] ? "1" : "0");
    }

    string res = "";
    for (int i = 0; i < cf.n3; ++i) {
        res += (out[i] ? "1" : "0");
    }
    cout << "party:" << party << ", input: " << input << endl;
    cout << "party:" << party << ", result: " << res << endl;

    delete io;
    return 0;
}