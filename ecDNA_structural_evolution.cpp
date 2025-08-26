
# include <iostream>
# include <fstream>
# include <sstream>
# include <stdlib.h>
# include <stdio.h>
# include <cmath>
# include <math.h>
# include <random>
# include <cstdlib>
# include <dirent.h>
# include <string>
# include <getopt.h>
# include <algorithm>
# include <initializer_list>
# include <vector>
# include <iterator>




using namespace std;



/*******************************************************************************/



// Define global variables
const int _maxsize = 1e5;

int seed, initial_copyNumber, Ntot, N_ecDNA_hot, iter, doubled_ecDNA_copyNumber, dying_cell_index, daughter_ecDNA_copyNumber1, daughter_ecDNA_copyNumber2, num_ecDNA_segments, daughter_1_current_index, daughter_2_current_index;
double t, r_birth_normalised, r_death_normalised, total_unnormalised_division_rate, total_unnormalised_death_rate, rand_double, cumulative_division_rate, selection_coeff, sigmoid_a, sigmoid_b, p_change_size, ec_length_sum;
bool verbose_flag, BIRTH, DEATH, added_to_occupancy_vector, allocated_to_daughter_1;
vector<int> daughter_1_ec, daughter_2_ec, daughter_1_ec_indices, mother_ec, mother_cell_indices;


/*******************************************************************************/



// Define a cell
class Cell
{
	public:
		vector<int> ecDNA;
		double division_rate;

		// Constructors for Cell object
		Cell()
		{
			vector<int> ec;
			set_ecDNA(ec);
			set_division_rate(ec);
		}

		Cell(vector<int> ec)
		{
			set_ecDNA(ec);
			set_division_rate(ec);
		}

		// Set() and get() methods
		
		void set_ecDNA(vector<int> ec)
		{
			this->ecDNA = ec;
		}

		void set_division_rate(vector<int> ec)
		{
			// Sigmoid selection function
			// Compute mean ecDNA length 
			if (ec.size() == 0) this->division_rate = 1.0;
			else
			{
				ec_length_sum = 0.0;
				for (int i = 0; i < ec.size(); ++i) ec_length_sum += (double)ec[i];
				this->division_rate = 1.0 + (1.0 - (((ec_length_sum/(double)ec.size()) - 1)/(double)(num_ecDNA_segments-1)))*(selection_coeff/(1.0+exp(-sigmoid_a * (ec.size() - sigmoid_b))));
			}

			// Constant selection function
			// if (ec.size() == 0) this->division_rate = 1.0;
			// else this->division_rate = 1.0 + selection_coeff;
		}
};








//-----------------------







// Define an element of the occupancy vector
class Occupancy
{
	public:
		int ecDNA;
		int multiplicity;

		// Constructors for occupancy object
		Occupancy(int i , int j)
		{
			set_ecDNA(i);
			set_multiplicity(j);
		}

		// Set() and get() methods
		
		void set_ecDNA(int n)
		{
			this->ecDNA = n;
		}

		void set_multiplicity(int n)
		{
			this->multiplicity = n;
		}


};





//-----------------------






// Select cell to divide, based on individual division rates
void cell_division(vector<Cell> &tissue , double *total_unnormalised_division_rate , int *Ntot , int *N_ecDNA_hot , int num_ecDNA_segments , double p_change_size , mt19937_64 *generator)
{

	cumulative_division_rate = 0.0;

	rand_double = drand48();

	for (int i = 0; i < *Ntot; ++i)
	{

		cumulative_division_rate += tissue[i].division_rate;

		if (rand_double <= cumulative_division_rate/(*total_unnormalised_division_rate))
		{


			// Every ecDNA is copied once
			doubled_ecDNA_copyNumber = tissue[i].ecDNA.size() * 2;


			// cout << "Mother cell ecDNA (pre-ecDNA replication): ";
			// for (int j = 0; j < tissue[i].ecDNA.size(); ++j)
			// {
			// 	cout << tissue[i].ecDNA[j] << ", ";
			// }
			// cout << endl;


			// mother_ec is a full, explicit list of ecDNA in the mother cell after ecDNA replication
			mother_ec.clear();
			mother_ec.insert(mother_ec.begin(), tissue[i].ecDNA.begin(), tissue[i].ecDNA.end() );
			mother_ec.insert(mother_ec.begin(), tissue[i].ecDNA.begin(), tissue[i].ecDNA.end() );


			// cout << "Mother cell ecDNA (post-ecDNA replication): ";
			// for (int j = 0; j < mother_ec.size(); ++j)
			// {
			// 	cout << mother_ec[j] << ", ";
			// }
			// cout << endl;




			// Evolve ecDNA sizes
			for (int j = 0; j < mother_ec.size(); ++j)
			{
				// Edge cases
				if (mother_ec[j] == 1)
				{
					if (drand48() <= p_change_size/2.0) mother_ec[j] += 1;
				}
				else if (mother_ec[j] == num_ecDNA_segments)
				{
					if (drand48() <= p_change_size/2.0) mother_ec[j] -= 1;
				}

				// non-edge cases
				else if (drand48() <= p_change_size)
				{
					if (drand48() < 0.5) mother_ec[j] += 1;
					else mother_ec[j] -= 1;
				}
			}



			// cout << "Mother cell ecDNA (after size changes): ";
			// for (int j = 0; j < mother_ec.size(); ++j)
			// {
			// 	cout << mother_ec[j] << ", ";
			// }
			// cout << endl;




			// Distribute ecDNA between two daughter cells according to binomial
			binomial_distribution<int> draw_new_ecDNA_copyNumber(doubled_ecDNA_copyNumber , 0.5);
			daughter_ecDNA_copyNumber1 = draw_new_ecDNA_copyNumber(*generator);
			daughter_ecDNA_copyNumber2 = doubled_ecDNA_copyNumber - daughter_ecDNA_copyNumber1;


			// Book-keeping
			*total_unnormalised_division_rate -= tissue[i].division_rate;


			// Allocate ecDNA to new daughter cells 
			daughter_1_ec.assign(daughter_ecDNA_copyNumber1 , 0);
			daughter_2_ec.assign(daughter_ecDNA_copyNumber2 , 0);


			// Define sequence of indices up to size of mother cell ecDNA vector 
			mother_cell_indices.assign(doubled_ecDNA_copyNumber , 0);
			for (int i = 0; i < doubled_ecDNA_copyNumber; ++i)
			{
				mother_cell_indices[i] = i;
			}


			// Randomly select n of these (n = daughter_ecDNA_copyNumber1)
			daughter_1_ec_indices.clear();
			sample(mother_cell_indices.begin(), mother_cell_indices.end(), back_inserter(daughter_1_ec_indices), daughter_ecDNA_copyNumber1, *generator);



			// Allocate ecDNA depending on if their index is contained within random selection vector 
			daughter_1_current_index = 0;
			daughter_2_current_index = 0;
			for (int i = 0; i < mother_ec.size(); ++i)
			{
				allocated_to_daughter_1 = false;
				for (int j = 0; j < daughter_1_ec_indices.size(); ++j)
				{
					if (daughter_1_ec_indices[j] == i)
					{
						daughter_1_ec[daughter_1_current_index] = mother_ec[i];
						daughter_1_current_index += 1;
						allocated_to_daughter_1 = true;
						break;
					}
				}

				// If index not in list of indices for daughter 1, ecDNA goes into daughter 2
				if (allocated_to_daughter_1 == false)
				{
					daughter_2_ec[daughter_2_current_index] = mother_ec[i];
					daughter_2_current_index += 1;
				}
			}


			// cout << "Daughter 1 cell ecDNA: ";
			// for (int j = 0; j < daughter_1_ec.size(); ++j)
			// {
			// 	cout << daughter_1_ec[j] << ", ";
			// }
			// cout << endl;



			// cout << "Daughter 2 cell ecDNA: ";
			// for (int j = 0; j < daughter_2_ec.size(); ++j)
			// {
			// 	cout << daughter_2_ec[j] << ", ";
			// }
			// cout << endl;


		


			// Create daughter cells
			tissue[i] = Cell(daughter_1_ec);
			tissue[*Ntot] = Cell(daughter_2_ec);



			// Book-keeping
			*total_unnormalised_division_rate += tissue[i].division_rate + tissue[*Ntot].division_rate;

			if ((daughter_ecDNA_copyNumber1 > 0) && (daughter_ecDNA_copyNumber2 > 0)) *N_ecDNA_hot += 1;
			*Ntot += 1;

			return;
		}
	}


	// If we get this far, there is a problem with the Gillespie rates. 
	cout << "Problem with Gillespie rates encountered when choosing cell to divide. Exiting..." << endl;
	exit(0);
	
}






//-----------------------






// Select cell to daeth, based on individual death rates
void cell_death(vector<Cell> &tissue , double *total_unnormalised_division_rate , int *Ntot , int *N_ecDNA_hot , mt19937_64 *generator)
{

	cumulative_division_rate = 0.0;

	rand_double = drand48();

	for (int i = 0; i < *Ntot; ++i)
	{

		cumulative_division_rate += tissue[i].division_rate;

		if (rand_double <= cumulative_division_rate/(*total_unnormalised_division_rate))
		{

			// Book-keeping
			if (tissue[i].ecDNA.size() > 0) *N_ecDNA_hot -= 1;
			*total_unnormalised_division_rate -= tissue[i].division_rate;

			// Cell dies. Replace with last cell in tissue vector so that all vector indices from 0 to Ntot-1 are occupied by cells
			tissue[i] = Cell(tissue[*Ntot-1].ecDNA);
			*Ntot -= 1;
			
			return;
		}
	}


	// If we get this far, there is a problem with the Gillespie rates. 
	cout << "Problem with Gillespie rates encountered when choosing cell to divide. Exiting..." << endl;
	exit(0);
	
}







//-----------------------








// Parse command line arguments (Flags and numerical arguments)
void parse_command_line_arguments(int argc, char** argv , bool *verbose_flag , int *seed , int *initial_copyNumber , double *selection_coeff , double *sigmoid_a , double *sigmoid_b , int *num_ecDNA_segments , double *p_change_size)
{
	int c;
	int option_index;
	char* arg_long = nullptr;
	int verbose = 0;

	static struct option long_options[] =
	{
		{"verbose", no_argument, &verbose, 1},
	}; 

	while ((c = getopt_long(argc, argv, "x:n:s:a:b:l:p:", long_options, &option_index)) != -1)
	switch (c)
	{
		case 0:
		{
			arg_long = optarg;
			break;
		}

		// Random seed
		case 'x':
			*seed = atoi(optarg);		
			break;


		// ecDNA copy number in initial cell
		case 'n':
			*initial_copyNumber = atoi(optarg);
			break;


		// ecDNA copy number in initial cell
		case 's':
			*selection_coeff = atof(optarg);
			break;


		// a parameter for sigmoid fitness function
		case 'a':
			*sigmoid_a = atof(optarg);
			break;


		// b parameter for sigmoid fitness function
		case 'b':
			*sigmoid_b = atof(optarg);
			break;


		// parameter defining number of segments into which each ecDNA is divided
		case 'l':
			*num_ecDNA_segments = atoi(optarg);
			break;


		// parameter defining number of segments into which each ecDNA is divided
		case 'p':
			*p_change_size = atof(optarg);
			break;


		case '?':
			if (optopt == 'c')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr,
				"Unknown option character `\\x%x'.\n",
				optopt);
		//return 1;
		default:
		abort ();
	}



	// Set boolean values for verbose flag
	if (verbose == 1) *verbose_flag = true;


	if (*initial_copyNumber < 0)
	{
		cout << "Initial copy number must be 0 or greater. Exiting." << endl;
		exit(0);
	}

	if (*num_ecDNA_segments <= 1)
	{
		cout << "Invalid ecDNA size parameter. Exiting..." << endl;
		exit(0);
	}

	if ((*p_change_size < 0) || (*p_change_size > 1))
	{
		cout << "Invalid ecDNA size evolution probability. Exiting..." << endl;
		exit(0);
	}
}






//-----------------------






// Set up tissue (i.e. array of cells)
vector<Cell> initialise_tissue(int _maxsize , double *total_unnormalised_division_rate , int *Ntot , int *N_ecDNA_hot , int initial_copyNumber , int num_ecDNA_segments)
{

	if (verbose_flag) cout << " " << endl;


	// Set up the vector of cells, called tissue
	vector<Cell> tissue(_maxsize); 
	if (verbose_flag) printf(" Initialising tissue... Done.\r");
	if (verbose_flag) cout << " " << endl;
		

	// Seed first tissue cell
	vector<int> initial_ec(initial_copyNumber , num_ecDNA_segments);
	tissue[0] = Cell(initial_ec);


	// Book-keeping
	*Ntot += 1;
	*N_ecDNA_hot += 1;
	*total_unnormalised_division_rate = tissue[0].division_rate;

	return tissue;
}






//-----------------------






// Compute normalised birth and death rates 
void compute_normalised_birth_and_death_rates(int Ntot , int N_ecDNA_hot , double total_unnormalised_division_rate , double *r_birth_normalised , double *r_death_normalised)
{

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Compute un-normalised reaction rates, constant death rate for all cells
	//total_unnormalised_death_rate = 0.5*Ntot;
	//total_unnormalised_death_rate = 0.0;

	// If you want cell death to be proportional to cell birth rate, compute the line below and invoke cell_death() method
	total_unnormalised_death_rate = 0.5*total_unnormalised_division_rate;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



	// Compute normalised reaction rates
	*r_birth_normalised = total_unnormalised_division_rate/(total_unnormalised_division_rate + total_unnormalised_death_rate);
	*r_death_normalised = total_unnormalised_death_rate/(total_unnormalised_division_rate + total_unnormalised_death_rate);
}






//-----------------------






// Choose next event in Gillespie algorithm
void choose_next_event(bool *BIRTH , bool *DEATH , double r_birth_normalised , double r_death_normalised)
{
	*BIRTH = false;
	*DEATH = false;

	rand_double = drand48();
	if (rand_double <= r_birth_normalised)
	{
		*BIRTH = true;
	}
	else if (rand_double <= r_birth_normalised + r_death_normalised)
	{
		*DEATH = true;
	}
	else
	{
		cout << "Problem with Gillespie rates..." << endl;
		exit(0);
	}
}








//-----------------------







// // Kill cell and remove from lattice
// void kill_cell(vector<Cell> &tissue , int *Ntot , int *N_ecDNA_hot , double *total_unnormalised_division_rate)
// {

// 	// Randomly select one cell to die (uniform probability)
// 	dying_cell_index = round((drand48() * (*Ntot)) - 0.5);

	
// 	// Book-keeping
// 	if (tissue[dying_cell_index].ecDNA > 0) *N_ecDNA_hot -= 1;
// 	*total_unnormalised_division_rate -= tissue[dying_cell_index].division_rate;

// 	// Cell dies. Replace with last cell in tissue vector so that all vector indices from 0 to Ntot-1 are occupied by cells
// 	tissue[dying_cell_index] = Cell(tissue[*Ntot-1].ecDNA);
// 	*Ntot -= 1;

// }











// /*******************************************************************************/















int main(int argc, char** argv)
{




	// Reset time and tissue size variables
	t = 0.0;
	Ntot = 0;
	N_ecDNA_hot = 0;
	total_unnormalised_division_rate = 0.0;




	//================== Parse command line arguments ====================//
	parse_command_line_arguments(argc , argv , &verbose_flag , &seed , &initial_copyNumber , &selection_coeff , &sigmoid_a , &sigmoid_b , &num_ecDNA_segments , &p_change_size);



	// Seed random number generator
	srand48(seed);
	mt19937_64 generator;
	generator.seed(seed);








	//================== Initialise tissue ====================//
	vector<Cell> tissue = initialise_tissue(_maxsize , &total_unnormalised_division_rate , &Ntot , &N_ecDNA_hot , initial_copyNumber , num_ecDNA_segments);








	//================== Simulate tissue growth ==================//
	iter = 0;

	do
	{



		// Timestamp for simulation optimisation purposes
		//clock_t iter_start_time, iter_end_time;

		//iter_start_time = clock();

		//if (iter == 20) exit(0);

		//cout << "---- iter = " << iter << " ---------- n_ecDNA_hot = " << N_ecDNA_hot << " ----------------------------" << endl;

		
		++iter;




		// Re-evaluate birth and death rates
		compute_normalised_birth_and_death_rates(Ntot , N_ecDNA_hot , total_unnormalised_division_rate , &r_birth_normalised , &r_death_normalised);




		// Choose birth or death based on normalised rates
		choose_next_event(&BIRTH , &DEATH , r_birth_normalised , r_death_normalised);




		// If division:
		if (BIRTH)
		{
			// Cell divides
			cell_division(tissue , &total_unnormalised_division_rate , &Ntot , &N_ecDNA_hot , num_ecDNA_segments , p_change_size , &generator);
		}


		// If death:
		if ((DEATH) && (Ntot > 1))
		{
			// Cell dies (if constant death rate used)
			//kill_cell(tissue , &Ntot , &N_ecDNA_hot , &total_unnormalised_division_rate);

			// Cell dies (if CN dependent death rate used)
			cell_death(tissue , &total_unnormalised_division_rate , &Ntot , &N_ecDNA_hot , &generator);
		}







		if (iter%1000 == 0)
		{			
			if (verbose_flag) cout << "Iteration #" << iter << " -- N = " << Ntot << " -- N_ecDNA_hot = " << N_ecDNA_hot << endl;
		
			// If all ecDNA lost from population, exit Gillespie loop and write output data
			if (N_ecDNA_hot == 0)
			{
				stringstream f;
				f.str("");
				f << "./RESULTS/Nmax=" << _maxsize << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_l=" << num_ecDNA_segments << "_p=" << p_change_size << "/seed=" << seed;
				DIR *dir = opendir(f.str().c_str());
				if(!dir)
				{
					f.str("");
					f << "mkdir -p ./RESULTS/Nmax=" << _maxsize << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_l=" << num_ecDNA_segments << "_p=" << p_change_size << "/seed=" << seed;
					system(f.str().c_str());
				}

				ofstream tissue_file;
				f.str("");
				f << "./RESULTS/Nmax=" << _maxsize << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_l=" << num_ecDNA_segments << "_p=" << p_change_size << "/seed=" << seed << "/tissue.csv";
				tissue_file.open(f.str().c_str());


				if (verbose_flag) cout << " " << endl;
				if (verbose_flag) cout << "Created output files..." << endl;


				// Loop over all cells 
				for (int i = 0; i < tissue.size(); ++i)
				{
					// Write -1 for cells with no ecDNA
					tissue_file << "-1" << endl;
				}


				exit(0);
			}
		}





	} while (Ntot < _maxsize);		// Exit once system has reached total size of _maxsize

	if (verbose_flag) cout << " " << endl;










	//================== Open data files & write final system data ==================//



	stringstream f;
	f.str("");
	f << "./RESULTS/Nmax=" << _maxsize << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_l=" << num_ecDNA_segments << "_p=" << p_change_size << "/seed=" << seed;
	DIR *dir = opendir(f.str().c_str());
	if(!dir)
	{
		f.str("");
		f << "mkdir -p ./RESULTS/Nmax=" << _maxsize << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_l=" << num_ecDNA_segments << "_p=" << p_change_size << "/seed=" << seed;
		system(f.str().c_str());
	}

	ofstream tissue_file;
	f.str("");
	f << "./RESULTS/Nmax=" << _maxsize << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_l=" << num_ecDNA_segments << "_p=" << p_change_size << "/seed=" << seed << "/tissue.csv";
	tissue_file.open(f.str().c_str());


	if (verbose_flag) cout << " " << endl;
	if (verbose_flag) cout << "Created output files..." << endl;



	
	// Loop over all cells 
	for (int i = 0; i < tissue.size(); ++i)
	{
		// Write -1 for cells with no ecDNA
		if (tissue[i].ecDNA.size() == 0) tissue_file << "-1" << endl;

		// Otherwise, write list of ecDNA sizes
		else
		{
			for (int j = 0; j < tissue[i].ecDNA.size(); ++j)
			{
				if (j == tissue[i].ecDNA.size() - 1) tissue_file << tissue[i].ecDNA[j] << endl;
				else tissue_file << tissue[i].ecDNA[j] << ",";
			}
		}
	}



	tissue_file.close();







	return 0;
}

















