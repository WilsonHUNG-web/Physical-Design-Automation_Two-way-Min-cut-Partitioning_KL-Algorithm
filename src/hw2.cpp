#include<iostream>
#include<vector>
#include<list>
#include<map>
#include<fstream>
#include<string>
#include<chrono>
#include<algorithm>
#include <omp.h>
using namespace std;


//all necessary global variables
int pmax = 0; //initialize Max. pin = 0
int total_size = 0; //initialize total cell size = 0, set in parse_cell
int cut_size = 0;  //initialize total cut_size = 0
bool debug = false;  //initialize debug mode OFF
int total_net_num = 0;  //initialize total net number = 0
int total_cell_num = 0;  //initialize total cell number = 0
float ratio_fac = 0; //initialize r for balance
int smax = 0;


class net {
public:
	net(int name) : n_name(name) { cinA = 0; cinB = 0; }
	void add_cell(int c_name) { this->list_cell.push_back(c_name); }
	int name() { return this->n_name; }
	int get_cinA() { return this->cinA; } //return cell number in A
	int get_cinB() { return this->cinB; }
	int get_cell_num_in_bucket(int parti) {
		if (parti == 0) return this->cinA;
		else if (parti == 1) return this->cinB;
		else { cout << "Wrong partitition in [get_cell_num_in_bucket]" << endl; return -1; }
	}


	void update_cinA(int c_num) { this->cinA = c_num; }
	void update_cinB(int c_num) { this->cinB = c_num; }
	void update_bucket_cellnum(int parti,int cell_num) { 
		if (parti == 0) this->cinA = cell_num;
		else if (parti == 1) this->cinB= cell_num;
	}

	void print_netinfo() {
		cout << "n" << n_name << " has ";
		unsigned int i = 0;
		while (i < list_cell.size()) { cout << "c" << list_cell[i] << " "; i++; }
		cout << endl;
	}
	vector<int> list_cell;

private:
	int cinA;
	int cinB;
	int n_name;
};

class cell {
public:
	cell() {

	}
	cell(int name) : c_name(name) {
		lock = false;
		gain = 0;
		c_partition = -1; //doesn't belong to any partition in the beginning
	}

	void add_net(int n_name) { this->list_net.push_back(n_name); }
	void lock_cell() { this->lock = true; }
	void free_cell() { this->lock = false; }
	void update_size(int size) { this->c_size = size; }
	void update_partition(int parti) { this->c_partition = parti; }
	void reverse_parition() {
		if (this->c_partition == 0) this->c_partition = 1;
		else this->c_partition = 0;
	}
	void zeroing_gain() { this->gain = 0; }
	void add_gain() { this->gain = this->gain +1 ; }
	void reduce_gain() { this->gain = this->gain - 1; }
	void revers_gain() { this->gain = -(this->gain); }

	int name() { return this->c_name; }
	bool get_lock() { return this->lock; }
	bool is_free() { return !(this->lock); }
	int get_size() { return this->c_size; }
	int get_gain() { return this->gain; }
	int get_partition() { return this->c_partition; }
	int get_pin() { return this->list_net.size(); }
	vector<int> list_net;
	list<cell*>::iterator ptr;

	void operator=(cell* c) {
		this->c_name = c->c_name;
		this->c_size = c->c_size;
		this->gain = c->gain;
		this->lock = c->lock;
		this->c_partition = c->c_partition;
		this->list_net = c->list_net;

	}
	int c_name;
private:

	int c_size;
	int gain;

	bool lock;
	int c_partition;
};

map<int, cell*> cell_dic;
map<int, net*> net_dic;



int c_num;
int n_num;

vector<net*> vector_n;
vector<cell*> vector_c;

chrono::high_resolution_clock::time_point time_record() {
	return chrono::high_resolution_clock::now();
}

double time_output(chrono::high_resolution_clock::time_point start_time, chrono::high_resolution_clock::time_point end_time) {
	return chrono::duration<double>(end_time - start_time).count();
}

chrono::high_resolution_clock::time_point begin_time = time_record();

void set_maxpin() {

	for (int i = 0; i < c_num; i++) {
		//cout << "c" << vector_c[i]->name() << "'s p = " << vector_c[i]->get_pin() << endl;
		int temp = vector_c[i]->get_pin();
		if (temp > pmax) pmax = temp;
	}

	if (debug == true) cout << "pmax = " << pmax << endl;
	return;
}
//create buckets
map<int, list<cell*>> buc_A, buc_B;
int size_A = 0;
int size_B = 0;
int A_num = 0;
int B_num = 0;

//count cell numbers in A and B for each net
void set_net_c_num_inAnB() {
	for (unsigned int i = 0; i < vector_n.size(); i++) //for each net, update its A B num
	{
		vector_n[i]->update_cinA(0); //INITIALIZED THE CIN_A
		vector_n[i]->update_cinB(0); //INITIALIZED THE CIN_B
		for (unsigned int j = 0; j < vector_n[i]->list_cell.size(); j++) {
			int parti = cell_dic[vector_n[i]->list_cell[j]]->get_partition();
			if (parti == 0) {
				vector_n[i]->update_cinA( 1+ vector_n[i]->get_cinA());
			}
			else if (parti == 1) {
				vector_n[i]->update_cinB(1 + vector_n[i]->get_cinB());
			}

			else
			{
				cout << "cell series# " << i << " not initialized well with its partition!" << endl;
			}
		}
	}
	if (debug == true) cout << "[end set_net_c_num]" << endl;
	return;
}

void free_all_cells() {
	for (int i = 0; i < c_num; i++) {
		vector_c[i]->free_cell();
	}
	return;
}

void initialize_gain() {
	if(debug==true) cout << "[start initialize_gain]" << endl;
	

	for (int i = 0; i < c_num; i++) {
		vector_c[i]->zeroing_gain();
		for (unsigned int j = 0; j < vector_c[i]->list_net.size(); j++) {
			int net_name = vector_c[i]->list_net[j];

			if (vector_c[i]->get_partition() == 0) //if cell is in A
			{
				if (net_dic[net_name]->get_cinA() == 1) {
					vector_c[i]->add_gain(); //gain++
				}
				if (net_dic[net_name]->get_cinB() == 0) {
					vector_c[i]->reduce_gain();  //gain--
				}
			}

			else if (vector_c[i]->get_partition() == 1) //if cell is in B
			{
				if (net_dic[net_name]->get_cinB() == 1) {
					vector_c[i]->add_gain(); //gain++
				}
				if (net_dic[net_name]->get_cinA() == 0) {
					vector_c[i]->reduce_gain(); //gain--
				}
			}
		}
	}
	
	if (debug == true) cout << "[end initialize_gain]" << endl;
	return;
}

void reset_buckets() {
	if (debug == true) cout << "[start set_buckets]" << endl;
	
	buc_A.clear();
	buc_B.clear();
	
	//create list of cell* for each gain in each bucket
	for(int i =0; i<2*pmax ; i++) {
		list<cell*> gain_A_list, gain_B_list;
		buc_A[i] = gain_A_list;
		buc_B[i] = gain_B_list;
		i++;
	}

	//push cell to bucket according to its partition: A or B & its gain

	for (int i = 0; i < c_num; i++) {
		int gain = vector_c[i]->get_gain() +pmax;
		int parti = vector_c[i]->get_partition();

		if (parti == 0) {
			buc_A[gain].push_back(vector_c[i]);
		}

		else if (parti == 1) {
			buc_B[gain].push_back(vector_c[i]);
		}

		else
		{
			cout << "cell series# " << i << " not initialized well with its partition!" << endl;
			return;
		}
	}
}


//ok (?)
void initialize_partition_by_cellorder() { //initializer index = 1
	int a_size = 0;//initialize size of A
	int b_size = 0;

//assign the fist half of cells to part A by the order of cell list
	int i = 0;
	while (a_size < (total_size / 2)) {
		vector_c[i]->update_partition(0); //cell partition changed to A
		a_size = a_size + vector_c[i]->get_size();
		i++;
	}
	A_num = i;
	size_A = a_size;

	//assign the rest of cells to part B
	while (i < c_num) {
		vector_c[i]->update_partition(1); //cell partition changed to B
		b_size = b_size + vector_c[i]->get_size();
		i++;
	}
	B_num = i - A_num;
	size_B = b_size;
}


//ok
void initialize_partition_by_netorder() { //initializer index = 2
	int a_size = 0;//initialize size of A
	int b_size = 0;
	//assign the fist half of cells to part A by the order of net list
	unsigned int i = 0;

	while (i < vector_n.size()) {
		
		int net_name = vector_n[i]->name();
		for (unsigned int j = 0; j < net_dic[net_name]->list_cell.size(); j++) {
			int cell_name = net_dic[net_name]->list_cell[j];
			if (cell_dic[cell_name]->get_partition() != (-1)) continue; //if the cell belongs to A or B, skip and do next j
			cell_dic[cell_name]->update_partition(0); //cell partition changed to A
			a_size = a_size + cell_dic[cell_name]->get_size();
			A_num++;
			if (a_size >= (total_size / 2)) break;
		}

		if (a_size >= (total_size / 2)) break;
		i++;
	}

	size_A = a_size;

	int size_left = total_size - a_size;

	while (i < vector_n.size()|| b_size == size_left){
		int net_name = vector_n[i]->name();
		for (unsigned int j = 0; j < net_dic[net_name]->list_cell.size(); j++) {
			int cell_name = net_dic[net_name]->list_cell[j];
			if (cell_dic[cell_name]->get_partition() != (-1)) continue; //if the cell belongs to A or B, stop action
			cell_dic[cell_name]->update_partition(1); //cell partition changed to B
			b_size = b_size + cell_dic[cell_name]->get_size();
			B_num++;
			if (b_size == size_left) break;
		}
		i++;
		if (b_size == size_left) break;
	}

	size_B = b_size;

}
bool cmp(const net* &front, const net* &after)
{
	return front->list_cell.size() < after->list_cell.size();
}

void initial_by_net_celllist() {
	int a_size = 0;//initialize size of A
	int b_size = 0;
	int n_num = vector_n.size();
	//int a_till=0;
	vector<net*> vector_n_sort;
	cout << "starts" << endl;
	for (int i = 0; i < n_num; i++)
	{
		vector_n_sort.push_back(vector_n[i]);
	}

	cout << "finished" << endl;
	cout << "starts" << endl;
	//smort(vector_n_sort.begin(), vector_n_sort.end(), cmp);
	sort(vector_n_sort.begin(), vector_n_sort.end(),
		[](const net* a, const net* b) {
		return a->list_cell.size() > b->list_cell.size();
	});

	cout << "finished" << endl;
	unsigned int i;
	while (i < vector_n_sort.size()) {

		int net_name = vector_n_sort[i]->name();

		for (unsigned int j = 0; j < net_dic[net_name]->list_cell.size(); j++) {
			int cell_name = net_dic[net_name]->list_cell[j];
			if (cell_dic[cell_name]->get_partition() != (-1)) continue; //if the cell belongs to A or B, skip and do next j
			cell_dic[cell_name]->update_partition(0); //cell partition changed to A
			a_size = a_size + cell_dic[cell_name]->get_size();
			A_num++;
			if (a_size >= (total_size / 2)) break;
		}

		if (a_size >= (total_size / 2)) break;
		i++;
	}

	size_A = a_size;

	int size_left = total_size - a_size;

	while (i < vector_n_sort.size() || b_size == size_left) {
		int net_name = vector_n_sort[i]->name();

		for (unsigned int j = 0; j < net_dic[net_name]->list_cell.size(); j++) {
			int cell_name = net_dic[net_name]->list_cell[j];
			if (cell_dic[cell_name]->get_partition() != (-1)) continue; //if the cell belongs to A or B, stop action
			cell_dic[cell_name]->update_partition(1); //cell partition changed to B
			b_size = b_size + cell_dic[cell_name]->get_size();
			B_num++;
			if (b_size == size_left) break;
		}
		i++;
		if (b_size == size_left) break;
	}

	size_B = b_size;

}


//ok
void initializer_partition() {
	if (debug == true) cout << "[initialize_partition_by_netorder]" << endl;
	initialize_partition_by_netorder();
	//initial_by_net_celllist();
	return;

}

void count_cell_num_inAB() {

	int tmp_a = 0;
	int tmp_b = 0;
	for (int i = 0; i < c_num; i++)
	{
		int parti = vector_c[i]->get_partition();
		if (parti == 0) tmp_a++;
		else if (parti == 1) tmp_b++;
	}

	A_num = tmp_a;  //overwrite A_num
	B_num = tmp_b;  //overwrite B_num

	//check if A+B=cell_num
	if ((A_num + B_num) != c_num) {
		cout << "A+B != total cell number, we have some bugs " << endl;
	}
	return;
}

//ok//parce nets by file path 
void parce_nets(string file_path) {
	ifstream input;
	input.open(file_path);
	if (!input.is_open()) {
		cout << "Failed to load *.NETS" << endl;
	}
	string token;
	while (input >> token) {
		int tmp_name = -1;
		while (input >> token) {
			if (token == "NET") continue;
			if (token == "{") break;
			token = token.erase(0, 1);
			tmp_name = stoi(token);
		}

		net* temp_n = new net(tmp_name);
		net_dic[tmp_name] = temp_n;
		vector_n.push_back(net_dic[tmp_name]);
		while (input >> token) {
			if (token == "}") break;
			token = token.erase(0, 1);
			int c_name = stoi(token);
			if (!cell_dic[c_name]) {
				cell* temp_c = new cell(c_name);
				cell_dic[c_name] = temp_c;

				temp_c->add_net(tmp_name);
			}
			else {
				cell_dic[c_name]->add_net(tmp_name);
			}
			temp_n->add_cell(c_name);
		}
	}
	if (debug == true) cout << "[Net all parced]" << endl;
	
	input.close();
	total_net_num = vector_n.size();
	return;
}

//return the cell with max. gain, searching from the bucket top(pmax)
int get_maxgain_cell(map<int, list<cell*>> bucket) {

	//if (debug == true) cout << "get_maxgain_cell" << endl;
	int cell_name = -1;

	for (int i = 2*pmax; i >= 0; i--) {
		if (bucket[i].empty()) { 
			continue; }
		for (list<cell*>::iterator iter = bucket[i].begin(); iter != bucket[i].end(); iter++) {
				cell_name = (*iter)->name(); //if 1st one is found, break for loop
				if (debug == true) cout << "max cell name:" << cell_name << endl;
				return cell_name;
			//}
		}
		
	}
	if (cell_name == -1) cout << "bucket is empty " << endl;
	
	return cell_name;
}

//parce cells by file path
void parce_cells(string file_path) {
	//if (debug == true) cout << "START parce_cells" << endl;
	int sum_size = 0;
	ifstream input;
	input.open(file_path);
	if (!input.is_open()) {
		cout << "Failed to load *.CELLS" << endl;
		return;
	}
	char dummy_char;
	int name;
	int size;
	while (input >> dummy_char >> name >> size)
	{
		//cout << "cell name:" << name <<" cell size:" << size << endl;

		cell_dic[name]->update_size(size); //update each cell size
		vector_c.push_back(cell_dic[name]); //put every cellname into a vector
		sum_size = sum_size + size; //summantion size
	}

	total_size = sum_size;
	c_num = vector_c.size();
	input.close();
	return;
}


//ok
void set_cutsize() {
	cut_size = 0;
	if (debug == true) cout << "[start set_cutsize]" << endl;
	for (int i = 0; i < total_net_num; i++) {
		if (vector_n[i]->get_cinA() != 0 && vector_n[i]->get_cinB() != 0) {
			cut_size += 1;
		}
	}
	return;
}

void output_file(string file_path) {
	ofstream output;

	output.open(file_path);
	set_cutsize();
	output << "cut_size " << cut_size << endl;
	output << "A " << A_num << endl;
	for (int i = 0; i < c_num; i++) {
		if (vector_c[i]->get_partition() == 0) {
			output << "c" << vector_c[i]->name() << endl;
		}
	}
	
	output << "B " << B_num << endl;
	for (int i = 0; i < c_num; i++) {
		if (vector_c[i]->get_partition() == 1) {
			output << "c" << vector_c[i]->name() << endl;
		}
	}

	output.close();
	return;
}

//erase base cell from bucket list
void erase_base_cell_from_bucket(int c_name) {

	int parti = cell_dic[c_name]->get_partition();  //shall we use & ?
	int gain = cell_dic[c_name]->get_gain() + pmax;

	if (parti == 1) {
		bool flag = false;
		for (list<cell*>::iterator iter= buc_A[gain].begin(); iter != buc_A[gain].end(); iter++) {
			if (*iter == cell_dic[c_name]) {
				buc_A[gain].erase(iter); //remove from A
				flag = true;
				break;
			}
		}
		if (!flag) cout << "[erase_base_cell_from_bucket] Cannot find base cell in bucket A" << endl;
	}
	else if (parti == 0)
	{
		bool flag = false;
		for (list<cell*>::iterator  iter = buc_B[gain].begin(); iter != buc_B[gain].end(); iter++) {
			if (*iter == cell_dic[c_name]) {
				buc_B[gain].erase(iter); //remove from B
				flag = true;
				break;
			}
		}
		if (!flag) cout << "[erase_from_bucket] Cannot find base cell in bucket B" << endl;
	}

	else { cout << "[erase_from_bucket] cell has wrong partition" << endl; }

	return;
}


void erase_from_bucket(int c_name) {
	int parti = cell_dic[c_name]->get_partition();  //shall we use & ?
	int gain = cell_dic[c_name]->get_gain() + pmax;

	if (parti == 0) {
		bool flag = false;
		for (list<cell*>::iterator iter = buc_A[gain].begin(); iter != buc_A[gain].end(); iter++) {
			if (*iter == cell_dic[c_name]) {
				buc_A[gain].erase(iter); //remove from A
				flag = true;
				break;
			}
		}
		if (!flag) cout << "[erase_from_bucket] Cannot find base cell in bucket A" << endl;
	}
	else if (parti == 1)
	{
		bool flag = false;
		for (list<cell*>::iterator iter = buc_B[gain].begin(); iter != buc_B[gain].end(); iter++) {
			if (*iter == cell_dic[c_name]) {
				buc_B[gain].erase(iter); //remove from B
				flag = true;
				break;
			}
		}
		if (!flag) cout << "[erase_from_bucket] Cannot find base cell in bucket B" << endl;
	}

	else { cout << "[erase_from_bucket] cell has wrong partition" << endl; }

	return;
}

int reverse_partition(int parti) {
	if (parti == 0) return 1;
	else if (parti == 1) return 0;
	else {
		cout << "reverse_partition has invalid input" << endl;
		return -1; 
	}
}

void place_to_bucket(int c_name) {
	int parti = cell_dic[c_name]->get_partition();

	if (parti == 0);


}

void print_num_locked() {
	int count = 0;
	for (int i = 0; i < c_num; i++) {
		if (vector_c[i]->get_lock()) count++;
	}
	cout << "locked cell#: " << count << " total cell#: " << c_num << endl;
}

//lock the moving cell and update all the gains, manage the buckets
void update_gain(int c_name) {

	cell_dic[c_name]->lock_cell();  //lock the base cell 
	int Front = cell_dic[c_name]->get_partition();  //Front block
	int To = reverse_partition(Front); //To block

	int net_count = cell_dic[c_name]->list_net.size();


	for (int i = 0; i < net_count; i++) { //for each net n on the base cell

		int net_name = cell_dic[c_name]->list_net[i];

		net* net_temp = net_dic[net_name]; //& removed
		int cell_list_size = net_temp->list_cell.size();

		//T(n)=0
		if(net_temp->get_cell_num_in_bucket(To) == 0){ //T(n)=0
				for (int j = 0; j < cell_list_size; j++) {  //each cell on net n
				int cell_name_on_net = net_temp->list_cell[j];
					if (cell_dic[cell_name_on_net]->is_free()) {   //if cell is free
						erase_from_bucket(cell_name_on_net);
						cell_dic[cell_name_on_net]->add_gain(); //gain++
						int new_gain = pmax + cell_dic[cell_name_on_net]->get_gain();
						if (cell_dic[cell_name_on_net]->get_partition() == 0) buc_A[new_gain].push_front(cell_dic[cell_name_on_net]);
						else buc_B[new_gain].push_front(cell_dic[cell_name_on_net]);

					}//move cell to new list by its new gain
				//}
			}
		}

		//T(n)=1
		else if (net_temp->get_cell_num_in_bucket(To) == 1) { //T(n)=1
			for (int j = 0; j < cell_list_size; j++) {  //each cell on net n
				int cell_name_on_net = net_temp->list_cell[j];
				if (cell_dic[cell_name_on_net]->get_partition() == To) {
					if (cell_dic[cell_name_on_net]->is_free()) {   //if only 1 cell in T not locked
						erase_from_bucket(cell_name_on_net);
						cell_dic[cell_name_on_net]->reduce_gain(); //gain--
						int new_gain = pmax + cell_dic[cell_name_on_net]->get_gain();
						if (cell_dic[cell_name_on_net]->get_partition() == 0) buc_A[new_gain].push_front(cell_dic[cell_name_on_net]);
						else buc_B[new_gain].push_front(cell_dic[cell_name_on_net]);
					}
					break; //if got one updated or not, leave for loop
				}
			}
		}


		cell_dic[c_name]->update_partition(To); //B to A or A to B
		//T(n)=T(n)+1
		int T_increment = net_temp->get_cell_num_in_bucket(To) + 1;
		net_temp->update_bucket_cellnum(To, T_increment);
		//F(n) = F(n) ¡V1
		int F_decrement = net_temp->get_cell_num_in_bucket(Front) - 1;
		net_temp->update_bucket_cellnum(Front, F_decrement);



		//F(n)=0
		if (net_temp->get_cell_num_in_bucket(Front) == 0) { //F(n)=0
			for (int j = 0; j < cell_list_size; j++) {  //each cell on net n
				int cell_name_on_net = net_temp->list_cell[j];
					if (cell_dic[cell_name_on_net]->is_free()) {   //if cell is free
						erase_from_bucket(cell_name_on_net);
						cell_dic[cell_name_on_net]->reduce_gain(); //gain--
						int new_gain = pmax + cell_dic[cell_name_on_net]->get_gain();
						if (cell_dic[cell_name_on_net]->get_partition() == 0) buc_A[new_gain].push_front(cell_dic[cell_name_on_net]);
						else buc_B[new_gain].push_front(cell_dic[cell_name_on_net]);
					}
				//}
			}
		}
		//F(n)=1
		else if (net_temp->get_cell_num_in_bucket(Front) == 1) { //F(n)=1
			for (int j = 0; j < cell_list_size; j++) {  //each cell on net n
				int cell_name_on_net = net_temp->list_cell[j];
				if (cell_dic[cell_name_on_net]->get_partition() == Front) { //only 1 cell in F
					if (cell_dic[cell_name_on_net]->is_free()) {   //if cell is free
						erase_from_bucket(cell_name_on_net);
						cell_dic[cell_name_on_net]->add_gain(); //gain++
						int new_gain = pmax + cell_dic[cell_name_on_net]->get_gain();
						if (cell_dic[cell_name_on_net]->get_partition() == 0) buc_A[new_gain].push_front(cell_dic[cell_name_on_net]);
						else buc_B[new_gain].push_front(cell_dic[cell_name_on_net]);
					}
					break; //if got one updated or not, leave for loop
				}
			}
		}

	}//end for loop


	erase_base_cell_from_bucket(c_name); //base cell leave bucket
	return;
}


int get_maxsize(int parti) {
	int smax = 0;
	for (int i = 0; i < c_num; i++) {
		if (vector_c[i]->get_partition() == parti) {
			int temp = vector_c[i]->get_size();
			if (temp > smax) smax = temp;
		}
	}
	return smax;
}

void set_smax() {
	for (int i = 0; i < c_num; i++) {
			int temp = vector_c[i]->get_size();
			if (temp > smax) smax = temp;
		}

return;
}


void set_ratio_factor() {
	ratio_fac = (float)size_A / (float)total_size;
}

//bool balance_check(int parti) {
bool balance_check(int temp_sizeA) {

	//int smax = get_maxsize(0);
	int W = total_size;  //total cell size

	float lower = ratio_fac * (float)W - float(smax);
	float upper = ratio_fac * (float)W + float(smax);

	if (temp_sizeA >= lower && temp_sizeA <= upper) { 
		if (debug == true) cout << "balance_check good" << endl;
		return true; }
	else {
		if (debug == true) cout << "balance_check failed" << endl;
		return false; }
	
}

int decide_action(int A_free, int B_free) {

	int c_name_A = get_maxgain_cell(buc_A); //return -1 if bucket is empty or all locked
	int c_name_B = get_maxgain_cell(buc_B);

	if (c_name_A == -1) { if (c_name_B != -1) { //cout << "(c_name_A == -1 || c_name_B != -1)" << endl; 
	return 2; } }
	if (c_name_B == -1) { if (c_name_A != -1) return 1; }

	if (c_name_B == -1 && c_name_A == -1) { cout << "(c_name_B == -1 && c_name_A == -1)" << endl; return 0; }

	if (B_free == 0) return 1;
	if (A_free == 0) { //cout << "(A_free == 0)" << endl;  
		return 2; }
	else if (cell_dic[c_name_A]->get_gain() >= cell_dic[c_name_B]->get_gain()) return 3;
	else return 4;
	return 0;
}


int partition_pass3() {

	int count = c_num;
	int A_free = A_num;
	int B_free = B_num;
	int key = 0;
	int part_sum=0;
	int max_part_sum=0;
	vector<cell*> moved_cells;
	bool flag = false;
	int iter=0;
	int iter_max_sum=0;
	int countdown = 10;

	while (count && !flag) {
		int cell_gmax_A=-1;
		int cell_gmax_B=-1;


		//decide_action by key
		if (B_free == 0 && A_free == 0) key = 0;
		else if (B_free == 0) key = 1;
		else if (A_free == 0) key = 2;
		else {
			cell_gmax_A = get_maxgain_cell(buc_A);  //cell name
			cell_gmax_B = get_maxgain_cell(buc_B);

			if (cell_gmax_A == -1 && cell_gmax_B == -1) key = 0;
			else if (cell_gmax_A != -1 && cell_gmax_B == -1) key = 1;
			else if (cell_gmax_A == -1 && cell_gmax_B != -1) key = 2;
			else {
				int gain_A = cell_dic[cell_gmax_A]->get_gain();
				int gain_B = cell_dic[cell_gmax_B]->get_gain();

				if (gain_A >= gain_B) key = 3; //move A
				else key = 4; //move B
			}
		}

		if (debug) cout << "  case key = " <<key << endl;


		if (key == 0) break; //end loop

		else if (key == 1) { //move A
			int temp_sizeA = size_A - cell_dic[cell_gmax_A]->get_size();
			int temp_sizeB = size_B + cell_dic[cell_gmax_A]->get_size();

			if (balance_check(temp_sizeA)) {
				size_A = temp_sizeA;
				size_B = temp_sizeB;
				part_sum = part_sum + cell_dic[cell_gmax_A]->get_gain();
				moved_cells.push_back(cell_dic[cell_gmax_A]);
				update_gain(cell_gmax_A);
			}
			else {
				flag = true;
				break;
			}

		}


		else if (key == 2) { //move B
			int temp_sizeA = size_A + cell_dic[cell_gmax_B]->get_size();
			int temp_sizeB = size_B - cell_dic[cell_gmax_B]->get_size();

			if (balance_check(temp_sizeA)) {
				//if (debug) { cout << "base cell: " << cell_gmax_B << endl; }
				size_A = temp_sizeA;
				size_B = temp_sizeB;
				part_sum = part_sum + cell_dic[cell_gmax_B]->get_gain();
				moved_cells.push_back(cell_dic[cell_gmax_B]);
				update_gain(cell_gmax_B);
			}
			else {
				flag = true;
				break;
			}
		}
		else if (key == 3) { //move A, if failed, move B
			int temp_sizeA = size_A - cell_dic[cell_gmax_A]->get_size();
			int temp_sizeB = size_B + cell_dic[cell_gmax_A]->get_size();

			if (balance_check(temp_sizeA)) {
				//if (debug) { cout << "base cell: " << cell_gmax_A << endl; }
				size_A = temp_sizeA;
				size_B = temp_sizeB;
				part_sum = part_sum + cell_dic[cell_gmax_A]->get_gain();
				moved_cells.push_back(cell_dic[cell_gmax_A]);
				update_gain(cell_gmax_A);
			}
			else { //A faild, move B
				temp_sizeA = size_A + cell_dic[cell_gmax_B]->get_size();
				temp_sizeB = size_B - cell_dic[cell_gmax_B]->get_size();
				if (balance_check(temp_sizeA)) {
					//if (debug) { cout << "base cell: " << cell_gmax_B << endl; }
					size_A = temp_sizeA;
					size_B = temp_sizeB;
					part_sum = part_sum + cell_dic[cell_gmax_B]->get_gain();
					moved_cells.push_back(cell_dic[cell_gmax_B]);
					update_gain(cell_gmax_B);
				}
				else { //both move A and B failed, end loop
					flag = true;
					break;
				}
			}
		}
		else { //key==4:  move B, if failed, move A
			int temp_sizeA = size_A + cell_dic[cell_gmax_B]->get_size();
			int temp_sizeB = size_B - cell_dic[cell_gmax_B]->get_size();

			if (balance_check(temp_sizeA)) {
				size_A = temp_sizeA;
				size_B = temp_sizeB;
				part_sum = part_sum + cell_dic[cell_gmax_B]->get_gain();
				moved_cells.push_back(cell_dic[cell_gmax_B]);
				update_gain(cell_gmax_B);
			}
			else { //B faild, move A
				temp_sizeA = size_A - cell_dic[cell_gmax_A]->get_size();
				temp_sizeB = size_B + cell_dic[cell_gmax_A]->get_size();
				if (balance_check(temp_sizeA)) {
					size_A = temp_sizeA;
					size_B = temp_sizeB;
					part_sum = part_sum + cell_dic[cell_gmax_A]->get_gain();
					moved_cells.push_back(cell_dic[cell_gmax_A]);
					update_gain(cell_gmax_A);
				}
				else { //both move B and A failed, end loop
					flag = true;
					break;
				}
			}
			if(debug) cout << "part_sum: " << part_sum << endl;
		}
		iter++;

		countdown--;
		if (part_sum >= max_part_sum) {
			max_part_sum = part_sum;
			iter_max_sum = iter;
			countdown++;
		}
		count--;
		if (countdown == 0) break;

		chrono::high_resolution_clock::time_point pass_time = time_record();
		if (time_output(begin_time, pass_time) >= 550) break;
	} //end of while loop

	//restore to max_part_sum_state by change back the partition
	int useless_iter = iter - iter_max_sum ;

	for (int i = 0; i < useless_iter; i++) {
		moved_cells.back()->reverse_parition();
		moved_cells.pop_back();
	}

	set_net_c_num_inAnB();  //update cell counts by A and B for each net coz we did restore action
	count_cell_num_inAB();  //recalculate global A_num and B_num coz we did restore action

	//if(debug)cout << "max_part_sum: " << max_part_sum << endl;
	return max_part_sum;
}



void partition() {
	if (debug == true) cout << "[start partition]" << endl;
	int count = 0;
	int i = 10;  //limit the max. pass #
	chrono::high_resolution_clock::time_point begin = time_record();
	chrono::high_resolution_clock::time_point end_time;

	while (i > 0) {
		int Gk;
		//if (debug == true) cout << "pass# " << count << "starts" << endl;
		Gk = partition_pass3(); //update Gk
		end_time = time_record();
		if (time_output(begin, end_time) >= 500) break;
		if (Gk <= 0)break;
			//cout << "[Gk <= 0] end partition" << endl;
			

		free_all_cells();
		initialize_gain();
		reset_buckets();
		
		i--;
	}
	


	if (debug == true) cout << "pass# " << count << "finished" << endl;
	
	
}

void debugMode(string arg) {
	if (arg == "-d") {
		debug = true;
		cout << "Debug Mode Turned ON!" << endl;
	}
}



int main(int argc, char *argv[]) {
	debugMode(argv[argc - 1]);
	//string path_prefix = "../testcases/";
	//if (!valid_command(argc)) return 0; //terminate if argc is not valid

	//if (debug == true) cout << "We have entered " << argc << " arguments"<<endl;

	string* file_path;
	file_path = new string[3];
	for (int i = 1; i <= 2; ++i) {
		string tmp_str(argv[i]);
		file_path[i] = tmp_str;

	}
	/*if (debug == true) {
		cout << "File path for nets: " << file_path[1] << endl;
		cout << "File path for cells: " << file_path[2] << endl;
	}*/

	//----parse start time----
	chrono::high_resolution_clock::time_point p_start = time_record();

	parce_nets(file_path[1]);
	parce_cells(file_path[2]);

	chrono::high_resolution_clock::time_point p_end = time_record();
	double p_time = time_output(p_start, p_end);
	//cout << "Parce time: " << p_time << " s" << endl;
	//----count parce time----


	//----COMPUTE start----
	chrono::high_resolution_clock::time_point comp_start = time_record();

	initializer_partition(); //assign cells to A or B partition by update_partition()
	set_maxpin();
	set_ratio_factor();
	set_smax();
	//construct_buckets(); //initializing bucket A + B with proper gain
	set_net_c_num_inAnB();
	initialize_gain();
	reset_buckets();
	partition();

	chrono::high_resolution_clock::time_point comp_end = time_record();
	double compute_time = time_output(comp_start, comp_end);
	//cout << "Compute time: " << compute_time << " s" << endl;
	//----count COMPUTE time----


	//----output start----
	chrono::high_resolution_clock::time_point out_start = time_record();

	output_file(argv[3]);

	chrono::high_resolution_clock::time_point out_end = time_record();
	double output_time = time_output(out_start, out_end);
	//cout << "output time: " << output_time << " s" << endl;
	//----count output time----



	cout << "IO time     :" << p_time + output_time << " s" << endl;
	cout << "Compute time: " << compute_time << " s" << endl;
	cout << "Total time  : " << p_time + output_time + compute_time << " s" << endl;

	return 0;
}