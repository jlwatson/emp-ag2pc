#ifndef EMP_AG2PC_METADATA_FILE_H
#define EMP_AG2PC_METADATA_FILE_H

#include "emp-tool/execution/protocol_execution.h"
#include "emp-tool/utils/block.h"
#include "emp-tool/circuits/bit.h"
#include <stdio.h>

namespace emp {
    class MetadataFile {
    public:
        int partition_id, num_input_wires, num_output_wires
        map<int, int> input_worker_map;
        map<int, int> output_worker_map;
        int tmp;
        MetadataFile(const char * file) {
            FILE *f = fopen(file, "r");
            tmp2 = fscanf(f, "%d%d%d\n", &partition_id, &num_input_wires, &num_output_wires);
            tmp2 = fscanf(f, "\n");
            for (int i = 0; i < num_input_wires; ++i) {
                int input_w, worker;
                tmp = fscanf(f, "%d%d", &input_w, &worker);
                input_worker_map.insert({input_w, worker})
            }
            tmp2 = fscanf(f, "\n");
            for (int i = 0; i < num_output_wires; ++i) {
                int out_w, worker;
                tmp = fscanf(f, "%d%d", &out_w, &worker);
                output_worker_map.insert({out_w, worker})
            }
            fclose(f);
        }
    };
}
#endif //EMP_AG2PC_METADATA_FILE_H
