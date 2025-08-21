
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




using namespace std;



/*******************************************************************************/



// Define global variables
const int _maxsize = 1e5;

#define PI 3.14159265


int seed, initial_copyNumber, Ntot, N_ecDNA_hot, iter, doubled_ecDNA_copyNumber, dividing_cell_index, dying_cell_index, daughter_ecDNA_copyNumber1, daughter_ecDNA_copyNumber2;
double t, r_birth_normalised, r_death_normalised, total_unnormalised_division_rate, total_unnormalised_death_rate, rand_double, cumulative_division_rate, selection_coeff, sigmoid_a, sigmoid_b;
bool verbose_flag, BIRTH, DEATH, added_to_occupancy_vector;




/*******************************************************************************/



// Define a cell
class Cell
{
	public:
		int ecDNA;
		double division_rate;

		// Constructors for Cell object
		Cell()
		{
			set_ecDNA(-1);
			set_division_rate(-1);
		}

		Cell(int n)
		{
			set_ecDNA(n);
			set_division_rate(n);
		}

		// Set() and get() methods
		
		void set_ecDNA(int n)
		{
			this->ecDNA = n;
		}

		void set_division_rate(int n)
		{
			this->division_rate = 1.0 + (selection_coeff/(1.0+exp(-sigmoid_a * (n - sigmoid_b))));
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
void cell_division(vector<Cell> &tissue , double *total_unnormalised_division_rate , int *Ntot , int *N_ecDNA_hot , mt19937_64 *generator)
{

	dividing_cell_index = -1;
	cumulative_division_rate = 0.0;

	rand_double = drand48();

	for (int i = 0; i < *Ntot; ++i)
	{

		cumulative_division_rate += tissue[i].division_rate;

		if (rand_double <= cumulative_division_rate/(*total_unnormalised_division_rate))
		{


			// Every ecDNA is copied once
			doubled_ecDNA_copyNumber = tissue[i].ecDNA * 2;


			// Distribute ecDNA between two daughter cells according to binomial
			binomial_distribution<int> draw_new_ecDNA_copyNumber(doubled_ecDNA_copyNumber , 0.5);
			daughter_ecDNA_copyNumber1 = draw_new_ecDNA_copyNumber(*generator);
			daughter_ecDNA_copyNumber2 = doubled_ecDNA_copyNumber - daughter_ecDNA_copyNumber1;


			// Book-keeping
			*total_unnormalised_division_rate -= tissue[i].division_rate;


			// Create daughter cell
			tissue[i] = Cell(daughter_ecDNA_copyNumber1);
			tissue[*Ntot] = Cell(daughter_ecDNA_copyNumber2);



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

	dividing_cell_index = -1;
	cumulative_division_rate = 0.0;

	rand_double = drand48();

	for (int i = 0; i < *Ntot; ++i)
	{

		cumulative_division_rate += tissue[i].division_rate;

		if (rand_double <= cumulative_division_rate/(*total_unnormalised_division_rate))
		{

			// Book-keeping
			if (tissue[i].ecDNA > 0) *N_ecDNA_hot -= 1;
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
void parse_command_line_arguments(int argc, char** argv , bool *verbose_flag , int *seed , int *initial_copyNumber, double *selection_coeff, double *sigmoid_a, double *sigmoid_b)
{
	int c;
	int option_index;
	char* arg_long = nullptr;
	int verbose = 0;

	static struct option long_options[] =
	{
		{"verbose", no_argument, &verbose, 1},
	}; 

	while ((c = getopt_long(argc, argv, "x:n:s:a:b:", long_options, &option_index)) != -1)
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
}






//-----------------------






// Set up tissue (i.e. array of cells)
vector<Cell> initialise_tissue(int _maxsize , double *total_unnormalised_division_rate , int *Ntot , int *N_ecDNA_hot , int initial_copyNumber)
{

	if (verbose_flag) cout << " " << endl;


	// Set up the vector of cells, called tissue
	vector<Cell> tissue(_maxsize); 
	if (verbose_flag) printf(" Initialising tissue... Done.\r");
	if (verbose_flag) cout << " " << endl;
		

	// Seed first tissue cell
	tissue[0] = Cell(initial_copyNumber);


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







// Kill cell and remove from lattice
void kill_cell(vector<Cell> &tissue , int *Ntot , int *N_ecDNA_hot , double *total_unnormalised_division_rate)
{

	// Randomly select one cell to die (uniform probability)
	dying_cell_index = round((drand48() * (*Ntot)) - 0.5);

	
	// Book-keeping
	if (tissue[dying_cell_index].ecDNA > 0) *N_ecDNA_hot -= 1;
	*total_unnormalised_division_rate -= tissue[dying_cell_index].division_rate;

	// Cell dies. Replace with last cell in tissue vector so that all vector indices from 0 to Ntot-1 are occupied by cells
	tissue[dying_cell_index] = Cell(tissue[*Ntot-1].ecDNA);
	*Ntot -= 1;

}











/*******************************************************************************/















int main(int argc, char** argv)
{




	// Reset time and tissue size variables
	t = 0.0;
	Ntot = 0;
	N_ecDNA_hot = 0;
	total_unnormalised_division_rate = 0.0;




	//================== Parse command line arguments ====================//
	parse_command_line_arguments(argc , argv , &verbose_flag , &seed , &initial_copyNumber , &selection_coeff, &sigmoid_a, &sigmoid_b);



	// Seed random number generator
	srand48(seed);
	mt19937_64 generator;
	generator.seed(seed);





	//================== Initialise tissue ====================//
	vector<Cell> tissue = initialise_tissue(_maxsize , &total_unnormalised_division_rate , &Ntot , &N_ecDNA_hot , initial_copyNumber);








	//================== Simulate tissue growth ==================//
	iter = 0;

	do
	{



		// Timestamp for simulation optimisation purposes
		//clock_t iter_start_time, iter_end_time;

		//iter_start_time = clock();

		//cout << iter << endl;
		
		++iter;




		// Re-evaluate birth and death rates
		compute_normalised_birth_and_death_rates(Ntot , N_ecDNA_hot , total_unnormalised_division_rate , &r_birth_normalised , &r_death_normalised);




		// Choose birth or death based on normalised rates
		choose_next_event(&BIRTH , &DEATH , r_birth_normalised , r_death_normalised);




		// If division:
		if (BIRTH)
		{
			// Cell divides
			cell_division(tissue , &total_unnormalised_division_rate , &Ntot , &N_ecDNA_hot , &generator);
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
				f << "./RESULTS/Nmax=" << _maxsize << "_initialCopyNumber=" << initial_copyNumber << "_s=" << selection_coeff << "_sigmoidA=" << sigmoid_a << "_sigmoidB=" << sigmoid_b << "/seed=" << seed;
				DIR *dir = opendir(f.str().c_str());
				if(!dir)
				{
					f.str("");
					f << "mkdir -p ./RESULTS/Nmax=" << _maxsize << "_initialCopyNumber=" << initial_copyNumber << "_s=" << selection_coeff << "_sigmoidA=" << sigmoid_a << "_sigmoidB=" << sigmoid_b << "/seed=" << seed;
					system(f.str().c_str());
				}

				ofstream tissue_file;
				f.str("");
				f << "./RESULTS/Nmax=" << _maxsize << "_initialCopyNumber=" << initial_copyNumber << "_s=" << selection_coeff << "_sigmoidA=" << sigmoid_a << "_sigmoidB=" << sigmoid_b << "/seed=" << seed << "/tissue.csv";
				tissue_file.open(f.str().c_str());


				if (verbose_flag) cout << " " << endl;
				if (verbose_flag) cout << "Created output files..." << endl;


				tissue_file << "0," << _maxsize << endl;

				exit(0);
			}
		}





	} while (Ntot < _maxsize);		// Exit once system has reached total size of _maxsize

	if (verbose_flag) cout << " " << endl;










	//================== Open data files & write final system data ==================//



	stringstream f;
	f.str("");
	f << "./RESULTS/Nmax=" << _maxsize << "_initialCopyNumber=" << initial_copyNumber << "_s=" << selection_coeff << "_sigmoidA=" << sigmoid_a << "_sigmoidB=" << sigmoid_b << "/seed=" << seed;
	DIR *dir = opendir(f.str().c_str());
	if(!dir)
	{
		f.str("");
		f << "mkdir -p ./RESULTS/Nmax=" << _maxsize << "_initialCopyNumber=" << initial_copyNumber << "_s=" << selection_coeff << "_sigmoidA=" << sigmoid_a << "_sigmoidB=" << sigmoid_b << "/seed=" << seed;
		system(f.str().c_str());
	}

	ofstream tissue_file;
	f.str("");
	f << "./RESULTS/Nmax=" << _maxsize << "_initialCopyNumber=" << initial_copyNumber << "_s=" << selection_coeff << "_sigmoidA=" << sigmoid_a << "_sigmoidB=" << sigmoid_b << "/seed=" << seed << "/tissue.csv";
	tissue_file.open(f.str().c_str());


	if (verbose_flag) cout << " " << endl;
	if (verbose_flag) cout << "Created output files..." << endl;



	
	// To reduce size of output file, compute occupancy vector 
	vector<Occupancy> occupancy_vector;

	// Loop over all cells in system and add to occupancy vector
	for (int i = 0; i < Ntot; i++)
	{

		added_to_occupancy_vector = false;

		// For each cell, check if occupancy object exists in vector
		for (int j = 0; j < occupancy_vector.size(); ++j)
		 {
		 	// If object exists in vector, add contribution from this cell
		 	if (occupancy_vector[j].ecDNA == tissue[i].ecDNA)
		 	{
		 		occupancy_vector[j].multiplicity += 1;
		 		added_to_occupancy_vector = true;
		 	}
		 } 	

		 // If object does not exist, create it and add contribution from this cell
		 if (added_to_occupancy_vector == false)
		 {
		 	occupancy_vector.push_back(Occupancy(tissue[i].ecDNA , 1));
		 }
	}




	// Write occupancy vector data to file
	for (int i = 0; i < occupancy_vector.size(); ++i)
	{
		tissue_file << occupancy_vector[i].ecDNA << "," << occupancy_vector[i].multiplicity << endl;
	}



	tissue_file.close();







	return 0;
}

















