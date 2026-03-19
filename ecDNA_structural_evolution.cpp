
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
int seed, initial_copyNumber, Ntot, N_ecDNA_hot, iter, doubled_ecDNA_copyNumber, dying_cell_index, daughter_ecDNA_copyNumber1, daughter_ecDNA_copyNumber2, num_ecDNA_segments, daughter_1_current_index, daughter_2_current_index, Nmax, occupancy, max_ecDNA_size;
double t, r_birth_normalised, r_death_normalised, total_unnormalised_division_rate, total_unnormalised_death_rate, rand_double, cumulative_division_rate , cumulative_death_rate, selection_coeff, sigmoid_a, sigmoid_b, p_change_size, ec_length_sum, x1, x2, x3, x4, ecDNA_size_multiplier, ecDNA_size_multiplier_factor;
bool verbose_flag, BIRTH, DEATH, added_to_occupancy_vector, allocated_to_daughter_1;
vector<int> daughter_1_ec, daughter_2_ec, daughter_1_ec_indices, mother_ec, mother_cell_indices, all_cell_indices, resampling_indices;


/*******************************************************************************/



// Define a cell
class Cell
{
	public:
		vector<int> ecDNA;
		double division_rate;
		double death_rate;

		// Constructors for Cell object
		Cell()
		{
			vector<int> ec;
			set_ecDNA(ec);
			set_division_rate(ec);
			set_death_rate(ec);
		}

		Cell(vector<int> ec)
		{
			set_ecDNA(ec);
			set_division_rate(ec);
			set_death_rate(ec);
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
			if (ec.size() == 0) 
			{
				this->division_rate = 1.0;
			}
			else if (ec.size() >= sigmoid_b)
			{
				this->division_rate = 1.0 + selection_coeff;
			}
			else
			{
				x1 = selection_coeff;
				x2 = 1.0 + ((sigmoid_b - (double)ec.size())/(sigmoid_b - sigmoid_a));
				x3 = ((double)ec.size() / sigmoid_b);
				x4 = (sigmoid_b / (sigmoid_b - sigmoid_a));

				this->division_rate = 1.0 + (x1 * x2 * (pow(x3 , x4)));
			}

		}


		void set_death_rate(vector<int> ec)
		{

			// Base death rate varies between 0.5->1.0 times the cell's birth rate, depending on total ecDNA burden
			ec_length_sum = 0.0;
			for (int i = 0; i < ec.size(); ++i) ec_length_sum += (double)ec[i];

			this->death_rate = 0.5 + (ecDNA_size_multiplier_factor/100.0*ec_length_sum);

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
void cell_division(vector<Cell> &tissue , double *total_unnormalised_division_rate , double *total_unnormalised_death_rate , int *Ntot , int *N_ecDNA_hot , int num_ecDNA_segments , double p_change_size , mt19937_64 *generator)
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



			// mother_ec is a full, explicit list of ecDNA in the mother cell after ecDNA replication
			mother_ec.clear();
			mother_ec.insert(mother_ec.begin(), tissue[i].ecDNA.begin(), tissue[i].ecDNA.end() );
			mother_ec.insert(mother_ec.begin(), tissue[i].ecDNA.begin(), tissue[i].ecDNA.end() );



			// Evolve ecDNA sizes
			for (int j = 0; j < mother_ec.size(); ++j)
			{
				// Edge cases
				if (mother_ec[j] == 1)
				{
					if (drand48() <= p_change_size/2.0) mother_ec[j] += 1;
				}

				// non-edge cases
				else if (drand48() <= p_change_size)
				{
					if (drand48() < 0.5) mother_ec[j] += 1;
					else mother_ec[j] -= 1;
				}
			}



			// Distribute ecDNA between two daughter cells according to binomial
			binomial_distribution<int> draw_new_ecDNA_copyNumber(doubled_ecDNA_copyNumber , 0.5);
			daughter_ecDNA_copyNumber1 = draw_new_ecDNA_copyNumber(*generator);
			daughter_ecDNA_copyNumber2 = doubled_ecDNA_copyNumber - daughter_ecDNA_copyNumber1;


			// Book-keeping
			*total_unnormalised_division_rate -= tissue[i].division_rate;
			*total_unnormalised_death_rate -= tissue[i].death_rate;


			// Allocate ecDNA to new daughter cells 
			daughter_1_ec.assign(daughter_ecDNA_copyNumber1 , 0);
			daughter_2_ec.assign(daughter_ecDNA_copyNumber2 , 0);


			// Define sequence of indices up to size of mother cell ecDNA vector 
			mother_cell_indices.assign(doubled_ecDNA_copyNumber , 0);
			for (int j = 0; j < doubled_ecDNA_copyNumber; ++j)
			{
				mother_cell_indices[j] = j;
			}


			// Randomly select n of these (n = daughter_ecDNA_copyNumber1)
			daughter_1_ec_indices.clear();
			sample(mother_cell_indices.begin(), mother_cell_indices.end(), back_inserter(daughter_1_ec_indices), daughter_ecDNA_copyNumber1, *generator);



			// Allocate ecDNA depending on if their index is contained within random selection vector 
			daughter_1_current_index = 0;
			daughter_2_current_index = 0;
			for (int k = 0; k < mother_ec.size(); ++k)
			{
				allocated_to_daughter_1 = false;
				for (int j = 0; j < daughter_1_ec_indices.size(); ++j)
				{
					if (daughter_1_ec_indices[j] == k)
					{
						daughter_1_ec[daughter_1_current_index] = mother_ec[k];
						daughter_1_current_index += 1;
						allocated_to_daughter_1 = true;
						break;
					}
				}

				// If index not in list of indices for daughter 1, ecDNA goes into daughter 2
				if (allocated_to_daughter_1 == false)
				{
					daughter_2_ec[daughter_2_current_index] = mother_ec[k];
					daughter_2_current_index += 1;
				}
			}



			// Create daughter cells
			tissue[i] = Cell(daughter_1_ec);
			tissue[*Ntot] = Cell(daughter_2_ec);



			// Book-keeping
			*total_unnormalised_division_rate += tissue[i].division_rate + tissue[*Ntot].division_rate;
			*total_unnormalised_death_rate += tissue[i].death_rate + tissue[*Ntot].death_rate;

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
void cell_death(vector<Cell> &tissue , double *total_unnormalised_division_rate , double *total_unnormalised_death_rate , int *Ntot , int *N_ecDNA_hot , mt19937_64 *generator)
{

	cumulative_death_rate = 0.0;

	rand_double = drand48();

	for (int i = 0; i < *Ntot; ++i)
	{

		cumulative_death_rate += tissue[i].death_rate;

		if (rand_double <= cumulative_death_rate/(*total_unnormalised_death_rate))
		{

			// Book-keeping
			if (tissue[i].ecDNA.size() > 0) *N_ecDNA_hot -= 1;
			*total_unnormalised_division_rate -= tissue[i].division_rate;
			*total_unnormalised_death_rate -= tissue[i].death_rate;

			// Cell dies. Replace with last cell in tissue vector so that all vector indices from 0 to Ntot-1 are occupied by cells
			tissue[i] = Cell(tissue[*Ntot-1].ecDNA);
			*Ntot -= 1;
			
			return;
		}
	}


	// If we get this far, there is a problem with the Gillespie rates. 
	cout << "Problem with Gillespie rates encountered when choosing cell to die. Exiting..." << endl;
	exit(0);
	
}







//-----------------------








// Parse command line arguments (Flags and numerical arguments)
void parse_command_line_arguments(int argc, char** argv , bool *verbose_flag , int *seed , int *Nmax , int *initial_copyNumber , double *selection_coeff , double *sigmoid_a , double *sigmoid_b , double *p_change_size , double *ecDNA_size_multiplier_factor)
{
	int c;
	int option_index;
	char* arg_long = nullptr;
	int verbose = 0;

	static struct option long_options[] =
	{
		{"verbose", no_argument, &verbose, 1},
	}; 

	while ((c = getopt_long(argc, argv, "x:n:k:s:a:b:p:c:", long_options, &option_index)) != -1)
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


		// Selection coefficient
		case 'n':
			*Nmax = atoi(optarg);
			break;


		// ecDNA copy number in initial cell
		case 'k':
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
		case 'p':
			*p_change_size = atof(optarg);
			break;


		// ecDNA burden penalty
		case 'c':
			*ecDNA_size_multiplier_factor = atof(optarg);		
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

	if ((*p_change_size < 0) || (*p_change_size > 1))
	{
		cout << "Invalid ecDNA size evolution probability. Exiting..." << endl;
		exit(0);
	}

	if (*sigmoid_a >= *sigmoid_b)
	{
		cout << "Must have a < b. Exiting..." << endl;
		exit(0);
	}
}






//-----------------------






// Set up tissue (i.e. array of cells)
vector<Cell> initialise_tissue(int Nmax , double *total_unnormalised_division_rate , double *total_unnormalised_death_rate , int *Ntot , int *N_ecDNA_hot , int initial_copyNumber , int num_ecDNA_segments , mt19937_64 *generator)
{

	if (verbose_flag) cout << " " << endl;


	// Set up the vector of cells, called tissue
	vector<Cell> tissue(Nmax); 
	if (verbose_flag) printf(" Initialising tissue... Done.\r");
	if (verbose_flag) cout << " " << endl;
		

	// Sample from distribution of initial ecDNA sizes to create list of ecDNA in initial cell
	normal_distribution draw_ecDNA_size{(double)(num_ecDNA_segments), 5.0};

	vector<int> initial_ec(initial_copyNumber , 0);
	for (int i = 0; i < initial_ec.size(); ++i)
	{
		initial_ec[i] = (int)(max(1.0 , draw_ecDNA_size(*generator)));
	}

	// Seed first tissue cell
	tissue[0] = Cell(initial_ec);


	// Book-keeping
	*Ntot += 1;
	*N_ecDNA_hot += 1;
	*total_unnormalised_division_rate = tissue[0].division_rate;
	*total_unnormalised_death_rate = tissue[0].death_rate;

	return tissue;
}






//-----------------------






// Compute normalised birth and death rates 
void compute_normalised_birth_and_death_rates(int Ntot , int N_ecDNA_hot , double total_unnormalised_division_rate , double total_unnormalised_death_rate , double *r_birth_normalised , double *r_death_normalised)
{

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













// /*******************************************************************************/















int main(int argc, char** argv)
{



	// Reset time and tissue size variables
	t = 0.0;
	Ntot = 0;
	N_ecDNA_hot = 0;
	total_unnormalised_division_rate = 0.0;
	total_unnormalised_death_rate = 0.0;
	num_ecDNA_segments = 10;

	double num_replatings = 3;
	int resampling_number = 10000;



	//================== Parse command line arguments ====================//
	parse_command_line_arguments(argc , argv , &verbose_flag , &seed , &Nmax , &initial_copyNumber , &selection_coeff , &sigmoid_a , &sigmoid_b , &p_change_size , &ecDNA_size_multiplier_factor);


	if (Nmax < resampling_number)
	{
		cout << "Resampling number must be larger than max population size. Exiting..." << endl;
		exit(0);
	}



	// Seed random number generator
	srand48(seed);
	mt19937_64 generator;
	generator.seed(seed);



	//================== Initialise tissue ====================//
	vector<Cell> tissue = initialise_tissue(Nmax , &total_unnormalised_division_rate , &total_unnormalised_death_rate , &Ntot , &N_ecDNA_hot , initial_copyNumber , num_ecDNA_segments , &generator);




	//================== Simulate tissue growth ==================//
	iter = 0;

	for (int repeat = 0; repeat < num_replatings; ++repeat)
	{

		//cout << repeat << endl;
		
		do
		{



			// Timestamp for simulation optimisation purposes
			//clock_t iter_start_time, iter_end_time;
			//iter_start_time = clock();
			//cout << "---- iter = " << iter << " ---------- n_ecDNA_hot = " << N_ecDNA_hot << " ----------------------------" << endl;

			
			++iter;




			// Re-evaluate birth and death rates
			compute_normalised_birth_and_death_rates(Ntot , N_ecDNA_hot , total_unnormalised_division_rate , total_unnormalised_death_rate , &r_birth_normalised , &r_death_normalised);




			// Choose birth or death based on normalised rates
			choose_next_event(&BIRTH , &DEATH , r_birth_normalised , r_death_normalised);




			// If division:
			if (BIRTH)
			{
				// Cell divides
				cell_division(tissue , &total_unnormalised_division_rate , &total_unnormalised_death_rate , &Ntot , &N_ecDNA_hot , num_ecDNA_segments , p_change_size , &generator);
			}


			// If death:
			if ((DEATH) && (Ntot > 1))
			{
				// Cell dies (if CN dependent death rate used)
				cell_death(tissue , &total_unnormalised_division_rate , &total_unnormalised_death_rate , &Ntot , &N_ecDNA_hot , &generator);
			}







			if (iter%1000 == 0)
			{			
				if (verbose_flag) cout << "Iteration #" << iter << " -- N = " << Ntot << " -- N_ecDNA_hot = " << N_ecDNA_hot << endl;
			
				// If all ecDNA lost from population, exit Gillespie loop and write output data
				if (N_ecDNA_hot != 0) continue;


				stringstream f;
				f.str("");
				f << "./RESULTS/Nmax=" << Nmax << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_C=" << ecDNA_size_multiplier_factor << "_p=" << p_change_size << "/seed=" << seed;
				DIR *dir = opendir(f.str().c_str());
				if(!dir)
				{
					f.str("");
					f << "mkdir -p ./RESULTS/Nmax=" << Nmax << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_C=" << ecDNA_size_multiplier_factor << "_p=" << p_change_size << "/seed=" << seed;
					system(f.str().c_str());
				}

				ofstream tissue_file;
				f.str("");
				f << "./RESULTS/Nmax=" << Nmax << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_C=" << ecDNA_size_multiplier_factor << "_p=" << p_change_size << "/seed=" << seed << "/tissue_resample" << repeat <<".csv";
				tissue_file.open(f.str().c_str());


				if (verbose_flag) cout << " " << endl;
				if (verbose_flag) cout << "Created output files..." << endl;


				// Loop over all cells 
				for (int i = 0; i < tissue.size(); ++i)
				{
					
					tissue_file << "0,0" << endl;
					
				}
				


				return 0;
			}





		} while (Ntot < Nmax);		// Exit once system has reached total size of Nmax





		//================== Open data files & write final system data ==================//
		stringstream f;
		f.str("");
		f << "./RESULTS/Nmax=" << Nmax << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_C=" << ecDNA_size_multiplier_factor << "_p=" << p_change_size << "/seed=" << seed;
		DIR *dir = opendir(f.str().c_str());
		if(!dir)
		{
			f.str("");
			f << "mkdir -p ./RESULTS/Nmax=" << Nmax << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_C=" << ecDNA_size_multiplier_factor << "_p=" << p_change_size << "/seed=" << seed;
			system(f.str().c_str());
		}

		ofstream tissue_file;
		f.str("");
		f << "./RESULTS/Nmax=" << Nmax << "_k=" << initial_copyNumber << "_s=" << selection_coeff << "_A=" << sigmoid_a << "_B=" << sigmoid_b << "_C=" << ecDNA_size_multiplier_factor << "_p=" << p_change_size << "/seed=" << seed << "/tissue_resample" << repeat << ".csv";
		tissue_file.open(f.str().c_str());


		if (verbose_flag) cout << " " << endl;
		if (verbose_flag) cout << "Created output files..." << endl;



		
		// Loop over all cells 
		for (int i = 0; i < tissue.size(); ++i)
		{
			// Find max ecDNA size in cell 
			max_ecDNA_size = 0;
			for (int j = 0; j < tissue[i].ecDNA.size(); ++j)
			{
				if (tissue[i].ecDNA[j] > max_ecDNA_size) max_ecDNA_size = tissue[i].ecDNA[j];
			}


			for (int size = 0; size <= max_ecDNA_size; ++size)                                                                                                                                                                                                       
			{
				occupancy = 0;
				for (int j = 0; j < tissue[i].ecDNA.size(); ++j)
				{
					if (tissue[i].ecDNA[j] == size) occupancy += 1;
				}

				if (occupancy == 0) continue;

				if (size < max_ecDNA_size)
				{
					tissue_file << size << "," << occupancy << ";";
				}
				else
				{
					tissue_file << size << "," << occupancy << endl;
				}
			}
		}
		tissue_file.close();





		// Sample small number of cells to repeat growth experiment
		all_cell_indices.assign(Nmax , 0);
		for (int i = 0; i < Nmax; ++i)
		{
			all_cell_indices[i] = i;
		}


		// Randomly select n of these (n = resampling_number)
		resampling_indices.clear();
		sample(all_cell_indices.begin(), all_cell_indices.end(), back_inserter(resampling_indices), resampling_number, generator);


		// Re-populate tissue vector and reset book keeping variables and regrow to Nmax cells
		vector<Cell> tissue_resample(resampling_number);
		Ntot = resampling_number;
		N_ecDNA_hot = 0;
		total_unnormalised_division_rate = 0.0;
		total_unnormalised_death_rate = 0.0;

		for (int i = 0; i < resampling_number; ++i)
		{
			tissue_resample[i] = tissue[resampling_indices[i]];
			
			if (tissue[resampling_indices[i]].ecDNA.size() > 0) N_ecDNA_hot += 1;
			total_unnormalised_division_rate += tissue[resampling_indices[i]].division_rate;
			total_unnormalised_death_rate+= tissue[resampling_indices[i]].death_rate;
		}

		for (int i = 0; i < resampling_number; ++i)
		{
			tissue[i] = tissue_resample[i];
		}




	}

	if (verbose_flag) cout << " " << endl;












	return 0;
}

















