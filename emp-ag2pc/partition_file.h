#ifndef PARTITION_FILE 
#define PARTITION_FILE

#include "emp-tool/execution/circuit_execution.h"
#include "emp-tool/execution/protocol_execution.h"
#include "emp-tool/utils/block.h"
#include "emp-tool/circuits/bit.h"
#include <stdio.h>
#include <vector>
#include <set>
#include <map>

namespace emp {
#define AND_GATE 0
#define XOR_GATE 1
#define NOT_GATE 2

class PartitionFile { 
public:
	int num_gate, num_wire, num_alice_inputs, num_bob_inputs, num_outputs;
    int num_partition_wire;
    int num_partition_gate = 0;
    int num_partition_and = 0;
    int num_partition_xor = 0;
    int num_partition_inv = 0;
	int *gates;
	block * wires;
    map<int, int> global_wire_idx_to_local;
    map<int, int> local_and_gate_idx_to_global;
	int tmp, tmp2;
	PartitionFile(const char * file) {
		FILE * f = fopen(file, "r");
		tmp2=fscanf(f, "%d%d\n", &num_gate, &num_wire);
		tmp2=fscanf(f, "%d%d%d\n", &num_alice_inputs, &num_bob_inputs, &num_outputs);
		tmp2=fscanf(f, "\n");

		char str[10];
        std::vector<int> gate_vec;
        std::set<int> wire_set;

        while (fscanf(f, "%d", &tmp) != EOF) {
            if (tmp == 2) {
                int input_wire1, input_wire2, output_wire, global_gate_idx;
                tmp2=fscanf(f, "%d%d%d%d%d%s", &tmp, &input_wire1, &input_wire2, &output_wire, &global_gate_idx, str);
                gate_vec.push_back(input_wire1);
                gate_vec.push_back(input_wire2);
                gate_vec.push_back(output_wire);
                wire_set.insert(input_wire1);
                wire_set.insert(input_wire2);
                wire_set.insert(output_wire);

				if (str[0] == 'A') {
                    local_and_gate_idx_to_global[num_partition_gate] = global_gate_idx;
                    gate_vec.push_back(AND_GATE);   
                    num_partition_and++;
                }
                else if (str[0] == 'X') {
                    // ignore global gate index
                    gate_vec.push_back(XOR_GATE);
                    num_partition_xor++;
                }
            } else if (tmp == 1) {
                int input_wire1, output_wire, global_gate_idx;
                tmp2 = fscanf(f, "%d%d%d%d%s", &tmp, &input_wire1, &output_wire, &global_gate_idx, str);
                gate_vec.push_back(input_wire1);
                gate_vec.push_back(0);
                gate_vec.push_back(output_wire);
                gate_vec.push_back(NOT_GATE);
                wire_set.insert(input_wire1);
                wire_set.insert(output_wire);
                local_gate_idx_to_global[num_partition_gate] = global_gate_idx;

                num_partition_inv++;
            }
            num_partition_gate++;
        }

        gates = &gate_vec[0];
        num_partition_wire = wire_set.size();
        wires = new block[num_partition_wire];
        counter = 0;
        for (auto wire : wire_set) {
            global_wire_idx_to_local[wire] = counter;
            counter++;
        }

		fclose(f);
	}

	PartitionFile(const PartitionFile& cf) {
		num_gate = cf.num_gate;
		num_wire = cf.num_wire;
        num_partition_gate = cf.num_partition_gate;
        num_partition_wire = cf.num_partition_wire;
		num_alice_inputs = cf.num_alice_inputs;
		num_bob_inputs = cf.num_bob_inputs;
		num_outputs = cf.num_outputs;
		gates = new int[num_partition_gate * 4];
		wires = new block[num_partition_wire];
		memcpy(gates, cf.gates, num_partition_gate*4*sizeof(int));
		memcpy(wires, cf.wires, num_partition_wire*sizeof(block));	
	}
	~PartitionFile(){
		delete[] gates;
		delete[] wires;
	}
	int table_size() const{
		return num_partition_gate*4;
	}

    // XXX: figure out if this is used anywhere
    /*
	void compute(block * out, block * in1, block * in2) {
		memcpy(wires, in1, n1*sizeof(block));
		memcpy(wires+n1, in2, n2*sizeof(block));
		for(int i = 0; i < num_gate; ++i) {
			if(gates[4*i+3] == AND_GATE) {
				wires[gates[4*i+2]] = CircuitExecution::circ_exec->and_gate(wires[gates[4*i]], wires[gates[4*i+1]]);
			}
			else if (gates[4*i+3] == XOR_GATE) {
				wires[gates[4*i+2]] = CircuitExecution::circ_exec->xor_gate(wires[gates[4*i]], wires[gates[4*i+1]]);
			}
			else  
				wires[gates[4*i+2]] = CircuitExecution::circ_exec->not_gate(wires[gates[4*i]]);
		}
		memcpy(out, &wires[num_wire-n3], n3*sizeof(block));
	}
    */
};
}
#endif// PARTITION_FILE
