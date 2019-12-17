#ifndef C2PCDist_H__
#define C2PCDist_H__
#include "fpre.h"
#include <emp-tool/emp-tool.h>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;
using std::flush;
using std::cout;
using std::endl;
//#define __debug
BOOST_IS_BITWISE_SERIALIZABLE(emp::block)


// NOTE(romilb): If serialization gives troubles, try the below code.
// BOOST_SERIALIZATION_SPLIT_FREE(emp::block)
// namespace boost { namespace serialization {
//        template<class Archive>
//        void save(Archive& ar, const emp::block& m, unsigned int) {
//            ar << m;
//        }
//        template<class Archive>
//        void load(Archive& ar, emp::block& m, unsigned int) {
//            ar >> m;
//        }
//    }} // namespace boost::serialization

namespace emp {
class C2PCDist { public:
    const static int SSP = 5;//5*8 in fact...
	const block MASK = makeBlock(0x0ULL, 0xFFFFFULL);
	Fpre* fpre = nullptr;
	block * mac = nullptr;
	block * key = nullptr;
	bool * value = nullptr;

	block * preprocess_mac = nullptr;
	block * preprocess_key = nullptr;
	bool* preprocess_value = nullptr;

	block * sigma_mac = nullptr;
	block * sigma_key = nullptr;
	bool * sigma_value = nullptr;

	block * labels = nullptr;

	bool * mask = nullptr;
	CircuitFile * cf;
	NetIO * io;
	int num_ands = 0;
	int party, total_pre;

    C2PCDist() {

    }

    C2PCDist(NetIO * io, int party, CircuitFile * cf) {
		this->party = party;
		this->io = io;
		this->cf = cf;
		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == AND_GATE)
				++num_ands;
		}
        cout << "Using distributed C2PC\n";
		cout << cf->n1<<" "<<cf->n2<<" "<<cf->n3<<" "<<num_ands<<endl<<flush;
		total_pre = cf->n1 + cf->n2 + num_ands;
		fpre = new Fpre(io, party, num_ands);

		key = new block[cf->num_wire];
		mac = new block[cf->num_wire];
		value = new bool[cf->num_wire];

		preprocess_mac = new block[total_pre];
		preprocess_key = new block[total_pre];
		preprocess_value = new bool[total_pre];

		//sigma values in the paper
		sigma_mac = new block[num_ands];
		sigma_key = new block[num_ands];
		sigma_value = new bool[num_ands];

		labels = new block[cf->num_wire];

		mask = new bool[cf->n1 + cf->n2];
	}
	~C2PCDist(){
		delete[] key;
		delete[] mac;
		delete[] value;
		delete[] mask;
		delete[] GT;
		delete[] GTK;
		delete[] GTM;
		delete[] GTv;

		delete[] preprocess_mac;
		delete[] preprocess_key;
		delete[] preprocess_value;

		delete[] sigma_mac;
		delete[] sigma_key;
		delete[] sigma_value;

		delete[] labels;
		delete fpre;
	}
	PRG prg;
	PRP prp;
	block (* GT)[4][2] = nullptr;
	block (* GTK)[4] = nullptr;
	block (* GTM)[4] = nullptr;
	bool (* GTv)[4] = nullptr;

	//not allocation
	block * ANDS_mac = nullptr;
	block * ANDS_key = nullptr;
	bool * ANDS_value = nullptr;

	void function_independent() {
		if(party == ALICE)
			prg.random_block(labels, cf->num_wire);

		fpre->refill();
		ANDS_mac = fpre->MAC;
		ANDS_key = fpre->KEY;
		ANDS_value = fpre->r;

		prg.random_bool(preprocess_value, total_pre);
		if(fpre->party == ALICE) {
			fpre->abit1[0]->send(preprocess_key, total_pre);
			fpre->io[0]->flush();
			fpre->abit2[0]->recv(preprocess_mac, preprocess_value, total_pre);
			fpre->io2[0]->flush();
		} else {
			fpre->abit1[0]->recv(preprocess_mac, preprocess_value, total_pre);
			fpre->io[0]->flush();
			fpre->abit2[0]->send(preprocess_key, total_pre);
			fpre->io2[0]->flush();
		}
		memcpy(key, preprocess_key, (cf->n1+cf->n2)*sizeof(block));
		memcpy(mac, preprocess_mac, (cf->n1+cf->n2)*sizeof(block));
		memcpy(value, preprocess_value, (cf->n1+cf->n2)*sizeof(bool));
	}

	void function_dependent() {
		int ands = cf->n1+cf->n2;
		bool * x1 = new bool[num_ands];
		bool * y1 = new bool[num_ands];
		bool * x2 = new bool[num_ands];
		bool * y2 = new bool[num_ands];

		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == AND_GATE) {
				key[cf->gates[4*i+2]] = preprocess_key[ands];
				mac[cf->gates[4*i+2]] = preprocess_mac[ands];
				value[cf->gates[4*i+2]] = preprocess_value[ands];
				++ands;
			}
		}

		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == XOR_GATE) {
				key[cf->gates[4*i+2]] = xorBlocks(key[cf->gates[4*i]], key[cf->gates[4*i+1]]);
				mac[cf->gates[4*i+2]] = xorBlocks(mac[cf->gates[4*i]], mac[cf->gates[4*i+1]]);
				value[cf->gates[4*i+2]] = logic_xor(value[cf->gates[4*i]], value[cf->gates[4*i+1]]);
				if(party == ALICE)
					labels[cf->gates[4*i+2]] = xorBlocks(labels[cf->gates[4*i]], labels[cf->gates[4*i+1]]);
			} else if (cf->gates[4*i+3] == NOT_GATE) {
				if(party == ALICE)
					labels[cf->gates[4*i+2]] = xorBlocks(labels[cf->gates[4*i]], fpre->Delta);
				value[cf->gates[4*i+2]] = value[cf->gates[4*i]];
				key[cf->gates[4*i+2]] = key[cf->gates[4*i]];
				mac[cf->gates[4*i+2]] = mac[cf->gates[4*i]];
			}
		}
		ands = 0;
		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == AND_GATE) {
				x1[ands] = logic_xor(value[cf->gates[4*i]], ANDS_value[3*ands]);
				y1[ands] = logic_xor(value[cf->gates[4*i+1]], ANDS_value[3*ands+1]);	
				ands++;
			}
		}
		if(party == ALICE) {
			send_bool(io, x1, num_ands);
			send_bool(io, y1, num_ands);
			recv_bool(io, x2, num_ands);
			recv_bool(io, y2, num_ands);
		} else {
			recv_bool(io, x2, num_ands);
			recv_bool(io, y2, num_ands);
			send_bool(io, x1, num_ands);
			send_bool(io, y1, num_ands);
		}
		for(int i = 0; i < num_ands; ++i) {
			x1[i] = logic_xor(x1[i], x2[i]); 
			y1[i] = logic_xor(y1[i], y2[i]); 
		}
		ands = 0;
		for(int i = 0; i < cf->num_gate; ++i) {
			if (cf->gates[4*i+3] == AND_GATE) {
				sigma_mac[ands] = ANDS_mac[3*ands+2];
				sigma_key[ands] = ANDS_key[3*ands+2];
				sigma_value[ands] = ANDS_value[3*ands+2];

				if(x1[ands]) {
					sigma_mac[ands] = xorBlocks(sigma_mac[ands], ANDS_mac[3*ands+1]);
					sigma_key[ands] = xorBlocks(sigma_key[ands], ANDS_key[3*ands+1]);
					sigma_value[ands] = logic_xor(sigma_value[ands], ANDS_value[3*ands+1]);
				}
				if(y1[ands]) {
					sigma_mac[ands] = xorBlocks(sigma_mac[ands], ANDS_mac[3*ands]);
					sigma_key[ands] = xorBlocks(sigma_key[ands], ANDS_key[3*ands]);
					sigma_value[ands] = logic_xor(sigma_value[ands], ANDS_value[3*ands]);
				}
				if(x1[ands] and y1[ands]) {
					if(party == ALICE)
						sigma_key[ands] = xorBlocks(sigma_key[ands], fpre->Delta);
					else
						sigma_value[ands] = not sigma_value[ands];
				}

#ifdef __debug
				block MM[] = {mac[cf->gates[4*i]], mac[cf->gates[4*i+1]], sigma_mac[ands]};
				block KK[] = {key[cf->gates[4*i]], key[cf->gates[4*i+1]], sigma_key[ands]};
				bool VV[] = {value[cf->gates[4*i]], value[cf->gates[4*i+1]], sigma_value[ands]};
				check(MM, KK, VV);
#endif
				ands++;
			}
		}//sigma_[] stores the and of input wires to each AND gates

		delete[] fpre->MAC;
		delete[] fpre->KEY;
		delete[] fpre->r;
		fpre->MAC = nullptr;
		fpre->KEY = nullptr;
		fpre->r = nullptr;
		GT = new block[num_ands][4][2];
		GTK = new block[num_ands][4];
		GTM = new block[num_ands][4];
		GTv = new bool[num_ands][4];
	
		ands = 0;
		block H[4][2];
		block K[4], M[4];
		bool r[4];
		for(int i = 0; i < cf->num_gate; ++i) {
			if(cf->gates[4*i+3] == AND_GATE) {
				r[0] = logic_xor(sigma_value[ands] , value[cf->gates[4*i+2]]);
				r[1] = logic_xor(r[0] , value[cf->gates[4*i]]);
				r[2] = logic_xor(r[0] , value[cf->gates[4*i+1]]);
				r[3] = logic_xor(r[1] , value[cf->gates[4*i+1]]);
				if(party == BOB) r[3] = not r[3];

				M[0] = xorBlocks(sigma_mac[ands], mac[cf->gates[4*i+2]]);
				M[1] = xorBlocks(M[0], mac[cf->gates[4*i]]);
				M[2] = xorBlocks(M[0], mac[cf->gates[4*i+1]]);
				M[3] = xorBlocks(M[1], mac[cf->gates[4*i+1]]);

				K[0] = xorBlocks(sigma_key[ands], key[cf->gates[4*i+2]]);
				K[1] = xorBlocks(K[0], key[cf->gates[4*i]]);
				K[2] = xorBlocks(K[0], key[cf->gates[4*i+1]]);
				K[3] = xorBlocks(K[1], key[cf->gates[4*i+1]]);
				if(party == ALICE) K[3] = xorBlocks(K[3], fpre->Delta);

				if(party == ALICE) {
					Hash(H, labels[cf->gates[4*i]], labels[cf->gates[4*i+1]], i, prp, fpre);
					for(int j = 0; j < 4; ++j) {
						H[j][0] = xorBlocks(H[j][0], M[j]);
						H[j][1] = xorBlocks(H[j][1], xorBlocks(K[j], labels[cf->gates[4*i+2]]));
						if(r[j]) 
							H[j][1] = xorBlocks(H[j][1], fpre->Delta);
#ifdef __debug
						check2(M[j], K[j], r[j]);
#endif
					}
					for(int j = 0; j < 4; ++j ) {
						send_partial_block<SSP>(io, &H[j][0], 1);
						io->send_block(&H[j][1], 1);
					}
				} else {
					memcpy(GTK[ands], K, sizeof(block)*4);
					memcpy(GTM[ands], M, sizeof(block)*4);
					memcpy(GTv[ands], r, sizeof(bool)*4);
#ifdef __debug
					for(int j = 0; j < 4; ++j)
						check2(M[j], K[j], r[j]);
#endif
					for(int j = 0; j < 4; ++j ) {
						recv_partial_block<SSP>(io, &GT[ands][j][0], 1);
						io->recv_block(&GT[ands][j][1], 1);
					}
				}
				++ands;
			}
		}
		delete[] x1;
		delete[] x2;
		delete[] y1;
		delete[] y2;

		block tmp;
		if(party == ALICE) {
			send_partial_block<SSP>(io, mac, cf->n1);
			for(int i = cf->n1; i < cf->n1+cf->n2; ++i) {
				recv_partial_block<SSP>(io, &tmp, 1);
				block ttt = xorBlocks(key[i], fpre->Delta);
				ttt =  _mm_and_si128(ttt, MASK);
				block mask_key = _mm_and_si128(key[i], MASK);
				tmp =  _mm_and_si128(tmp, MASK);
				if(block_cmp(&tmp, &mask_key, 1))
					mask[i] = false;
				else if(block_cmp(&tmp, &ttt, 1))
					mask[i] = true;
				else cout <<"no match! ALICE\t"<<i<<endl;
			}
		} else {
			for(int i = 0; i < cf->n1; ++i) {
				recv_partial_block<SSP>(io, &tmp, 1);
				block ttt = xorBlocks(key[i], fpre->Delta);
				ttt =  _mm_and_si128(ttt, MASK);
				tmp =  _mm_and_si128(tmp, MASK);
				block mask_key = _mm_and_si128(key[i], MASK);
				if(block_cmp(&tmp, &mask_key, 1)) {
					mask[i] = false;
				} else if(block_cmp(&tmp, &ttt, 1)) {
					mask[i] = true;
				}
				else cout <<"no match! BOB\t"<<i<<endl;
			}

			send_partial_block<SSP>(io, mac+cf->n1, cf->n2);
		}
	}

	void online (bool * input, bool * output) {
	    uint8_t * mask_input = new uint8_t[cf->num_wire];
		memset(mask_input, 0, cf->num_wire);
		block tmp;
#ifdef __debug
		for(int i = 0; i < cf->n1+cf->n2; ++i)
			check2(mac[i], key[i], value[i]);
#endif
		if(party == ALICE) {
			for(int i = cf->n1; i < cf->n1+cf->n2; ++i) {
				mask_input[i] = logic_xor(input[i - cf->n1], value[i]);
				mask_input[i] = logic_xor(mask_input[i], mask[i]);
			}
			io->recv_data(mask_input, cf->n1);
			io->send_data(mask_input+cf->n1, cf->n2);
			for(int i = 0; i < cf->n1 + cf->n2; ++i) {
				tmp = labels[i];
				if(mask_input[i]) tmp = xorBlocks(tmp, fpre->Delta);
				io->send_block(&tmp, 1);
			}
			//send output mask data
			send_partial_block<SSP>(io, mac+cf->num_wire - cf->n3, cf->n3);
		} else {
			for(int i = 0; i < cf->n1; ++i) {
				mask_input[i] = logic_xor(input[i], value[i]);
				mask_input[i] = logic_xor(mask_input[i], mask[i]);
			}
			io->send_data(mask_input, cf->n1);
			io->recv_data(mask_input+cf->n1, cf->n2);
			io->recv_block(labels, cf->n1 + cf->n2);
		}
		int ands = 0;

		// Evaluate circuit
		if(party == BOB) {
			for(int i = 0; i < cf->num_gate; ++i) {
				if (cf->gates[4*i+3] == XOR_GATE) {
					labels[cf->gates[4*i+2]] = xorBlocks(labels[cf->gates[4*i]], labels[cf->gates[4*i+1]]);
					mask_input[cf->gates[4*i+2]] = logic_xor(mask_input[cf->gates[4*i]], mask_input[cf->gates[4*i+1]]);
				} else if (cf->gates[4*i+3] == AND_GATE) {
					int index = 2*mask_input[cf->gates[4*i]] + mask_input[cf->gates[4*i+1]];
					block H[2];
					Hash(H, labels[cf->gates[4*i]], labels[cf->gates[4*i+1]], i, index, prp);
					GT[ands][index][0] = xorBlocks(GT[ands][index][0], H[0]);
					GT[ands][index][1] = xorBlocks(GT[ands][index][1], H[1]);

					block ttt = xorBlocks(GTK[ands][index], fpre->Delta);
					ttt =  _mm_and_si128(ttt, MASK);
					GTK[ands][index] =  _mm_and_si128(GTK[ands][index], MASK);
					GT[ands][index][0] =  _mm_and_si128(GT[ands][index][0], MASK);

					if(block_cmp(&GT[ands][index][0], &GTK[ands][index], 1))
						mask_input[cf->gates[4*i+2]] = false;
					else if(block_cmp(&GT[ands][index][0], &ttt, 1))
						mask_input[cf->gates[4*i+2]] = true;
					else 	cout <<ands <<"no match GT!"<<endl;
					mask_input[cf->gates[4*i+2]] = logic_xor(mask_input[cf->gates[4*i+2]], GTv[ands][index]);

					labels[cf->gates[4*i+2]] = xorBlocks(GT[ands][index][1], GTM[ands][index]);
					ands++;
				} else {
					mask_input[cf->gates[4*i+2]] = not mask_input[cf->gates[4*i]];	
					labels[cf->gates[4*i+2]] = labels[cf->gates[4*i]];
				}
			}
		}

		// Unmask
		if (party == BOB) {
			bool * o = new bool[cf->n3];
			for(int i = 0; i < cf->n3; ++i) {
				block tmp;
				recv_partial_block<SSP>(io, &tmp, 1);
				tmp =  _mm_and_si128(tmp, MASK);

				block ttt = xorBlocks(key[cf->num_wire - cf-> n3 + i], fpre->Delta);
				ttt =  _mm_and_si128(ttt, MASK);
				key[cf->num_wire - cf-> n3 + i] =  _mm_and_si128(key[cf->num_wire - cf-> n3 + i], MASK);

				if(block_cmp(&tmp, &key[cf->num_wire - cf-> n3 + i], 1))
					o[i] = false;
				else if(block_cmp(&tmp, &ttt, 1))
					o[i] = true;
				else 	cout <<"no match output label!"<<endl;
			}
			for(int i = 0; i < cf->n3; ++i) {
				output[i] = logic_xor(o[i], mask_input[cf->num_wire - cf->n3 + i]);
				output[i] = logic_xor(output[i], value[cf->num_wire - cf->n3 + i]);
			}
			delete[] o;
		}
		delete[] mask_input;
	}

	void check(block * MAC, block * KEY, bool * r, int length = 1) {
		if (party == ALICE) {
			io->send_data(r, length*3);
			io->send_block(&fpre->Delta, 1);
			io->send_block(KEY, length*3);
			block DD;io->recv_block(&DD, 1);

			for(int i = 0; i < length*3; ++i) {
				block tmp;io->recv_block(&tmp, 1);
				if(r[i]) tmp = xorBlocks(tmp, DD);
				if (!block_cmp(&tmp, &MAC[i], 1))
					cout <<i<<"\tWRONG ABIT!"<<endl<<flush;
			}

		} else {
			bool tmp[3];
			for(int i = 0; i < length; ++i) {
				io->recv_data(tmp, 3);
				bool res = (logic_xor(tmp[0], r[3*i] )) and (logic_xor(tmp[1], r[3*i+1]));
				if(res != logic_xor(tmp[2], r[3*i+2]) ) {
					cout <<i<<"\tWRONG!"<<endl<<flush;
				}
			}
			block DD;io->recv_block(&DD, 1);

			for(int i = 0; i < length*3; ++i) {
				block tmp;io->recv_block(&tmp, 1);
				if(r[i]) tmp = xorBlocks(tmp, DD);
				if (!block_cmp(&tmp, &MAC[i], 1))
					cout <<i<<"\tWRONG ABIT!"<<endl<<flush;
			}
			io->send_block(&fpre->Delta, 1);
			io->send_block(KEY, length*3);
		}
		io->flush();
	}

	void check2(block & MAC, block & KEY, bool  &r) {
		if (party == ALICE) {
			io->send_block(&fpre->Delta, 1);
			io->send_block(&KEY, 1);
			block DD;io->recv_block(&DD, 1);
			for(int i = 0; i < 1; ++i) {
				block tmp;io->recv_block(&tmp, 1);
				if(r) tmp = xorBlocks(tmp, DD);
				if (!block_cmp(&tmp, &MAC, 1))
					cout <<i<<"\tWRONG ABIT!2"<<endl<<flush;
			}
		} else {
			block DD;io->recv_block(&DD, 1);
			for(int i = 0; i < 1; ++i) {
				block tmp;io->recv_block(&tmp, 1);
				if(r) tmp = xorBlocks(tmp, DD);
				if (!block_cmp(&tmp, &MAC, 1))
					cout <<i<<"\tWRONG ABIT!2"<<endl<<flush;
			}
			io->send_block(&fpre->Delta, 1);
			io->send_block(&KEY, 1);
		}
		io->flush();
	}

	static void Hash(block H[4][2], const block & a, const block & b, uint64_t i, PRP prp, Fpre* fpre) {
		block A[2], B[2];
		A[0] = a; A[1] = xorBlocks(a, fpre->Delta);
		B[0] = b; B[1] = xorBlocks(b, fpre->Delta);
		A[0] = double_block(A[0]);
		A[1] = double_block(A[1]);
		B[0] = double_block(double_block(B[0]));
		B[1] = double_block(double_block(B[1]));

		H[0][1] = H[0][0] = xorBlocks(A[0], B[0]);
		H[1][1] = H[1][0] = xorBlocks(A[0], B[1]);
		H[2][1] = H[2][0] = xorBlocks(A[1], B[0]);
		H[3][1] = H[3][0] = xorBlocks(A[1], B[1]);
		for(uint64_t j = 0; j < 4; ++j) {
			H[j][0] = xorBlocks(H[j][0], _mm_set_epi64x(4*i+j, 0ULL));
			H[j][1] = xorBlocks(H[j][1], _mm_set_epi64x(4*i+j, 1ULL));
		}
		prp.permute_block((block *)H, 8);
	}

	static void Hash(block H[2], block a, block b, uint64_t i, uint64_t row, PRP prp) {
		a = double_block(a);
		b = double_block(double_block(b));
		H[0] = H[1] = xorBlocks(a, b);
		H[0] = xorBlocks(H[0], _mm_set_epi64x(4*i+row, 0ULL));
		H[1] = xorBlocks(H[1], _mm_set_epi64x(4*i+row, 1ULL));
		prp.permute_block((block *)H, 2);
	}

	static bool logic_xor(bool a, bool b) {
		return a!= b;
	}
	string tostring(bool a) {
		if(a) return "T";
		else return "F";
	}

    void online2 (bool * input, bool * output) {
        uint8_t * mask_input = new uint8_t[cf->num_wire];
        memset(mask_input, 0, cf->num_wire);
        #ifdef __debug
                for(int i = 0; i < cf->n1+cf->n2; ++i)
                    check2(mac[i], key[i], value[i]);
        #endif
        send_recv_masks(input, mask_input);

        // Evaluate circuit
        if(party == BOB) {
            bob_parallel_evaluate(input, mask_input, cf->num_gate, cf->gates,
                                  party, labels,
                                  GT, GTK, GTM, GTv, fpre->Delta);

            // Unmask
            bob_unmask_output(output, mask_input);
        }
        delete[] mask_input;
    }

    void send_recv_masks(const bool *input, uint8_t *mask_input) const {
        block tmp;
        if(party == ALICE) {
            for(int i = cf->n1; i < cf->n1 + cf->n2; ++i) {
                mask_input[i] = logic_xor(input[i - cf->n1], value[i]);
                mask_input[i] = logic_xor(mask_input[i], mask[i]);
            }
            io->recv_data(mask_input, cf->n1);
            io->send_data(mask_input + cf->n1, cf->n2);
            for(int i = 0; i < cf->n1 + cf->n2; ++i) {
                tmp = labels[i];
                if(mask_input[i]) tmp = xorBlocks(tmp, fpre->Delta);
                io->send_block(&tmp, 1);
            }
            //send output mask data
            send_partial_block<SSP>(io, mac + cf->num_wire - cf->n3, cf->n3);
        } else {
            for(int i = 0; i < cf->n1; ++i) {
                mask_input[i] = logic_xor(input[i], value[i]);
                mask_input[i] = logic_xor(mask_input[i], mask[i]);
            }
            io->send_data(mask_input, cf->n1);
            io->recv_data(mask_input + cf->n1, cf->n2);
            io->recv_block(labels, cf->n1 + cf->n2);
        }
    }

    void bob_unmask_output(bool *output, const uint8_t *mask_input) const {
        if (party == BOB) {
            bool * o = new bool[cf->n3];
            for(int i = 0; i < cf->n3; ++i) {
                block tmp;
                recv_partial_block<SSP>(io, &tmp, 1);
                tmp =  _mm_and_si128(tmp, MASK);

                block ttt = xorBlocks(key[cf->num_wire - cf-> n3 + i], fpre->Delta);
                ttt =  _mm_and_si128(ttt, MASK);
                key[cf->num_wire - cf-> n3 + i] =  _mm_and_si128(key[cf->num_wire - cf-> n3 + i], MASK);

                if(block_cmp(&tmp, &key[cf->num_wire - cf-> n3 + i], 1))
                    o[i] = false;
                else if(block_cmp(&tmp, &ttt, 1))
                    o[i] = true;
                else 	cout <<"no match output label!"<<endl;
            }
            for(int i = 0; i < cf->n3; ++i) {
                output[i] = logic_xor(o[i], mask_input[cf->num_wire - cf->n3 + i]);
                output[i] = logic_xor(output[i], value[cf->num_wire - cf->n3 + i]);
            }
            delete[] o;
        }
    }

    static void evaluate_gate(int i, int &ands, block fpreDelta, PartitionFile *partition_file, block *labels, uint8_t *mask_input, PRP &prp,
                                  block (* GT)[4][2], block (* GTK)[4], block (* GTM)[4], bool (* GTv)[4]) {

        const block MASK = makeBlock(0x0ULL, 0xFFFFFULL);

        int *gates = partition_file->gates;
        if (gates[4 * i + 3] == XOR_GATE) {
            labels[gates[4 * i + 2]] = xorBlocks(labels[gates[4 * i]], labels[gates[4 * i + 1]]);
            mask_input[gates[4 * i + 2]] = logic_xor(mask_input[gates[4 * i]],
                                                     mask_input[gates[4 * i + 1]]);
        } else if (gates[4 * i + 3] == AND_GATE) {
            int index = 2 * mask_input[gates[4 * i]] + mask_input[gates[4 * i + 1]];
            block H[2];
            Hash(H, labels[gates[4 * i]], labels[gates[4 * i + 1]], i, index, prp);
            GT[ands][index][0] = xorBlocks(GT[ands][index][0], H[0]);
            GT[ands][index][1] = xorBlocks(GT[ands][index][1], H[1]);

            block ttt = xorBlocks(GTK[ands][index], fpreDelta);
            ttt = _mm_and_si128(ttt, MASK);
            GTK[ands][index] = _mm_and_si128(GTK[ands][index], MASK);
            GT[ands][index][0] = _mm_and_si128(GT[ands][index][0], MASK);

            if (block_cmp(&GT[ands][index][0], &GTK[ands][index], 1))
                mask_input[gates[4 * i + 2]] = false;
            else if (block_cmp(&GT[ands][index][0], &ttt, 1))
                mask_input[gates[4 * i + 2]] = true;
            else cout << ands << "no match GT!" << endl;
            mask_input[gates[4 * i + 2]] = logic_xor(mask_input[gates[4 * i + 2]], GTv[ands][index]);

            labels[gates[4 * i + 2]] = xorBlocks(GT[ands][index][1], GTM[ands][index]);
            ands++;
        } else {
            mask_input[gates[4 * i + 2]] = not mask_input[gates[4 * i]];
            labels[gates[4 * i + 2]] = labels[gates[4 * i]];
        }
    }

    static void bob_parallel_evaluate(bool * input, uint8_t * partition_mask_input, PartitionFile *partition_file, MetadataFile *meta_file, int party,
                                  block *partition_labels, 
                                  block (* GT)[4][2], block (* GTK)[4], block (* GTM)[4], bool (* GTv)[4],
                                  block fpreDelta, mpi::environment &env, mpi::communicator &world) {

	    // updates mask_input and labels
        if(party == BOB) {
            PRP prp;
            /*
            #ifdef __debug
                for(int i = 0; i < cf->n1+cf->n2; ++i)
                    check2(mac[i], key[i], value[i]);
            #endif
            */

            int ands = 0;
            int num_gates_processed = 0;

            std::map<int, mpi::request> input_requests;
            std::map<int, std::vector<int>> pending_wires;

            for (int i = 0; i < partition_file->num_partition_gate; i++) {
                int input_wire1 = partition_file->gates[4*i];
                if (input_wire1 >= partition_file->num_alice_inputs + partition_file->num_bob_inputs) {
                    pending_wires[i].push_back(input_wire1);
                    // send request to worker
                    if (meta_file->input_worker_map.contains(input_wire1) && !input_requests.contains(input_wire1)) { 
                        input_requests[input_wire1] = world.irecv(meta_file->input_worker_map[input_wire1], MPI_INPUT_BIT_REQUEST, &input[input_wire1]);
                    }
                }

                int type = partition_file->gates[4*i + 3];
                if (type == AND_GATE) {
                    int input_wire2 = partition_file->gates[4*i + 1];
                    if (input_wire2 >= partition_file->num_alice_inputs + partition_file->num_bob_inputs) {
                        pending_wires[i].push_back(input_wire2);
                        // send request to worker
                        if (meta_file->input_worker_map.contains(input_wire2) && !input_requests.contains(input_wire2)) { 
                            input_requests[input_wire2] = world.irecv(meta_file->input_worker_map[input_wire2], MPI_INPUT_BIT_REQUEST, &input[input_wire2]);
                        }
                    }
                }
            }

            while (num_gates_processed < partition_file->num_partition_gate) {
                //cout<<"Processing gate " << i << endl;
                //fflush(stdout);
                for (auto &wire_request_pair : input_requests) {
                    mpi::request r = wire_request_pair.second;
                    if (r.test()) {
                        for (auto &gate_pending_bits_pair : pending_wires) {
                            auto pos = std::find(gate_pending_bits_pair.second.begin(), gate_pending_bits_pair.second.end(), wire_request_pair.first);
                            if (pos != gate_pending_bits_pair.second.end()) {
                                gate_pending_bits_pair.second.erase(pos);
                            }
                        }
                    }
                }

                for (auto &gate_pending_bits_pair : pending_wires) {
                    if (gate_pending_bits_pair.second.empty()) { // evaluate the gate

                        evaluate_gate(gate_pending_bits_pair.first, ands, fpreDelta, partition_file, partition_labels, partition_mask_input, prp,
                                GT, GTK, GTM, GTv);

                        num_gates_processed++;

                        // remove output bits from pending map
                        int output_wire = partition_file->gates[4*gate_pending_bits_pair.first + 2];
                        for (auto &gate_pending_bits_pair2 : pending_wires) {
                            auto pos = std::find(gate_pending_bits_pair2.second.begin(), gate_pending_bits_pair2.second.end(), output_wire);
                            if (pos != gate_pending_bits_pair2.second.end()) {
                                gate_pending_bits_pair2.second.erase(pos);
                            }
                        }

                        // send output bit to other workers
                        if (meta_file->output_worker_map.contains(output_wire)) {
                            world.isend(meta_file->output_worker_map[output_wire], MPI_INPUT_BIT_REQUEST, &input[output_wire]);
                        }
                    }
                }
            }
        }
    }
};

    class C2PCDist_state { public:
        block fpreDelta;
        int cf_num_gates;
        int cf_num_wires;
        int num_ands;
        std::vector<int> gates;
        int party;
        std::vector<block> labels;

        // n-d arrays flattened to 1d
        std::vector<block> GTK;
        std::vector<block> GTM;
        std::vector<bool> GTv;
        std::vector<block> GT;

        friend class boost::serialization::access;

        C2PCDist_state(){

        }

        C2PCDist_state(C2PCDist c2pc) {
            fpreDelta = c2pc.fpre->Delta;
            cf_num_gates = c2pc.cf->num_gate;
            cf_num_wires = c2pc.cf->num_gate;
            num_ands = c2pc.num_ands;
            gates = std::vector<int>(c2pc.cf->gates, c2pc.cf->gates + cf_num_gates*4);
            party = c2pc.party;
            labels = std::vector<block>(c2pc.labels, c2pc.labels + cf_num_wires);
            GTK = std::vector<block>(*c2pc.GTK, *c2pc.GTK + 4*num_ands);
            GTM = std::vector<block>(*c2pc.GTM, *c2pc.GTM + 4*num_ands);
            GTv = std::vector<bool>(*c2pc.GTv, *c2pc.GTv + 4*num_ands);
            GT = std::vector<block>(**c2pc.GT, **c2pc.GT + 4*2*num_ands);
        }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & fpreDelta;
            ar & cf_num_gates;
            ar & cf_num_wires;
            ar & num_ands;
            ar & gates;
            ar & party;
            ar & labels;
            ar & GTK;
            ar & GTM;
            ar & GTv;
            ar & GT;
        }
    };
}
#endif// C2PCDist_H__
