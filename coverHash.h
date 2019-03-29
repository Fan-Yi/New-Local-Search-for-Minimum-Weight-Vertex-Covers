#ifndef COVER_HASH
#define COVER_HASH

#include <bitset>

#include "constants.h"

//hash(C)= (\sum_{v_i \in C}2^i) mod p, where p = PRIME_NUM

class CoverHash{
private:
	LL* two_power_mod_p;
	//* hash_entries;
	LL curr_hash_entry;
	bitset<PRIME_NUM> hash_entries;
	//int hash_bound;
public:
	CoverHash(LL v_num)
	{
		two_power_mod_p = new LL[v_num + 2];
		two_power_mod_p[0] = 1;
		for(LL v = 1; v <= v_num; v++)
		{
			two_power_mod_p[v] = (two_power_mod_p[v-1] * 2) % PRIME_NUM;
		}
///////////////////////////////////////////////
		curr_hash_entry = 0;
		//hash_bound = 1;
	}
	~CoverHash()
	{
		delete[] two_power_mod_p;
	}
/*
	int get_hash_entry()
	{
		return curr_hash_entry;
	}
*/
	void forget_all_visits()
	{
		hash_entries.reset();
	}

	void mark_hash_entry()
	{
		hash_entries.set(curr_hash_entry);
	}
/////////////////////////////////////////////////

	bool curr_hash_entry_visited()
	{
		return hash_entries.test(curr_hash_entry);
	}

	void update_hash_wrt_add(LL v)
	{
		curr_hash_entry = (curr_hash_entry + two_power_mod_p[v]) % PRIME_NUM;
	}

	void update_hash_wrt_remove(LL v)
	{
		curr_hash_entry = (curr_hash_entry + PRIME_NUM - two_power_mod_p[v]) % PRIME_NUM;
	}

};

#endif
