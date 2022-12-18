#include<vector>
#include<bitset>
#include<random>
#include<iostream>
#include<ctime>
#include<chrono>
#include<time.h>
#include<fstream>
#include<iomanip>
#define run(X,Y) for(int i = 1; i <= X; ++i) Y;


using namespace std;
using namespace std::chrono;

typedef double ld;
typedef long long int ll;
typedef pair<ld, ld> lld;

#define  BIT_LEN 3000

const ld eps = 1e-5;
const ld eps_dif = 1e-10;
const ll precision = 5;
const ld pi = 3.141592653;

///HI LO and Dimensions get modified here
ld minim = numeric_limits<ld>::max();
bitset<BIT_LEN> bits;
mt19937 random_generator(time(NULL));

ld inline log_pow(ld a, ll b) {
	ld r = 1;
	while (b) {
		if (b & 1)
			r *= a;
		a *= a;
		b /= 2;
	}
	return r;
}

inline void random_bit_string(bitset<BIT_LEN>& x, ll total_length) {
	for (ll i = 0; i < total_length; ++i) {
		x[i] = random_generator() % 2;
	}
}

double rand01()
{
	return static_cast<double>(random_generator()) / (RAND_MAX + 1.0);
}


inline void
bits_to_double(bitset<BIT_LEN>& bits, vector<ld>& double_vect, ld LO, ld HI, ll bit_length, ll total_length) {
	ll poz = 0;
	for (ll i = 0; i < total_length; i += bit_length) {
		ll xi = 0LL;
		for (ll j = i; j < i + bit_length; ++j) {
			xi *= 2;
			xi += bits[j];
		}
		ld value = (ld)xi / (log_pow(2.000, bit_length) - 1.0000);
		double_vect[poz++] = LO + value * (HI - LO);
		//        cout << double_vect[poz - 1] << ' ';
	}
	//    cout << '\n';
}


ld evaluate(bitset<BIT_LEN>& candidate, ld(*func)(vector<ld>& x, ll d), ll dim, ld LO, ld HI, ll bit_length,
	ll total_length) {
	vector<ld> vct(dim, 0.0000000);
	bits_to_double(candidate, vct, LO, HI, bit_length, total_length);
	return func(vct, dim);
}

inline ld Rastrigin(vector<ld>& x, ll d) {
	ld fx = 10.00 * d;
	for (ll i = 0; i < d; ++i) {
		fx += x[i] * x[i] - cos(2 * pi * x[i]) * 10;
	}

	return fx;
}


inline ld De_Jongs(vector<ld>& x, ll d) {
	ld fx = 0;
	for (ll i = 0; i < d; ++i) {
		fx += x[i] * x[i];
	}

	return fx;
}

inline ld Michalewicz(vector<ld>& x, ll d) {
	ld fx = 0;
	for (ll i = 0; i < d; ++i) {
		fx -= sin(x[i]) * log_pow(sin((i + 1) * x[i] * x[i] / pi), 20);
	}

	return fx;
}

inline ld Schwefels(vector<ld>& x, ll d) {
	ld fx = d * 418.9829;
	for (ll i = 0; i < d; ++i) {
		fx -= x[i] * sin(sqrt(abs(x[i])));
	}
	return fx;
}

inline ld t1p_function(vector <ld>& x, ll d)
{
	return x[0] * x[0] * x[0] - 60 * x[0] * x[0] + 900 * x[0] + 100;
}

bitset<BIT_LEN>
best_improvement(bitset<BIT_LEN>& x, ld(*func)(vector<ld>& x, ll d), ll dim, vector<ld>& doubles, ld LO, ld HI,
	ll bit_length, ll total_length) {
	ld fx = evaluate(x, func, dim, LO, HI, bit_length, total_length);
	ld fx_best = fx;
	bitset<BIT_LEN> best = x;
	for (ll i = 0; i < total_length; ++i) {
		x[i] = !x[i];
		ld fx_act = evaluate(x, func, dim, LO, HI, bit_length, total_length);

		if (fx_act < fx_best) {
			fx_best = fx_act;
			best = x;
		}
		x[i] = !x[i];
	}
	return best;
}


bitset<BIT_LEN>
first_improvement(bitset<BIT_LEN>& x, ld(*func)(vector<ld>& x, ll d), ll dim, vector<ld>& doubles, ld LO, ld HI,
	ll bit_length, ll total_length) {
	ld fx = evaluate(x, func, dim, LO, HI, bit_length, total_length);
	ld fx_best = fx;
	bitset<BIT_LEN> best = x;

	vector<int> bit_nr;
	for (ll i = 0; i < total_length; ++i)
		bit_nr.push_back(i);
	random_shuffle(bit_nr.begin(), bit_nr.end());
	for (ll idx = 0; idx < total_length; ++idx) {
		int i = bit_nr[idx];
		x[i] = !x[i];
		ld fx_act = evaluate(x, func, dim, LO, HI, bit_length, total_length);
		if (fx_act < fx_best) {
			fx_best = fx_act;
			best = x;
			x[i] = !x[i];
			return best;
		}
		x[i] = !x[i];
	}
	return best;
}


bitset<BIT_LEN>
worst_improvement(bitset<BIT_LEN>& x, ld(*func)(vector<ld>& x, ll d), ll dim, vector<ld>& doubles, ld LO, ld HI,
	ll bit_length, ll total_length) {
	ld fx = evaluate(x, func, dim, LO, HI, bit_length, total_length);
	ld fx_best = numeric_limits<ld>::min();
	bitset<BIT_LEN> best = x;
	for (ll i = 0; i < total_length; ++i) {
		x[i] = !x[i];
		ld fx_act = evaluate(x, func, dim, LO, HI, bit_length, total_length);
		if (fx_act > fx_best && fx_act < fx) {
			fx_best = fx_act;
			best = x;
		}
		x[i] = !x[i];
	}
	if (fx_best == numeric_limits<ld>::min())
		return x;
	return best;
}

void Hill_Climber(ld(*func)(vector<ld>& x, ll d),
	bitset<BIT_LEN>(*improve)(bitset<BIT_LEN>& x, ld(*func)(vector<ld>& x, ll d), ll dim,
		vector<ld>& doubles, ld LO, ld HI, ll bit_length, ll total_length), ll dim,
	ll ITERATIONS, ld LO, ld HI, const char* file_name) {
	ld best = numeric_limits<ld>::max();
	ll bit_length = ceil(log2(pow(10, precision) * (HI - LO)));
	ll total_length = bit_length * dim;
	vector<ld> doubles(dim), sol(dim);
	ofstream f(file_name, std::ios_base::app);

	f << setprecision(10) << fixed;

	auto start = steady_clock::now();
	for (ll it = 0; it < ITERATIONS; ++it) {
		//        cout << it << '\n';
		bool local = false;
		bitset<BIT_LEN> candidate;
		random_bit_string(candidate, total_length);
		ld fx_candidate = evaluate(candidate, func, dim, LO, HI, bit_length, total_length);
		do {
			bitset<BIT_LEN> neighbour = improve(candidate, func, dim, doubles, LO, HI, bit_length, total_length);
			ld fx_neighbour = evaluate(neighbour, func, dim, LO, HI, bit_length, total_length);
			if (fx_neighbour < fx_candidate) {
				candidate = neighbour;
				fx_candidate = fx_neighbour;
			}
			else local = true;
		} while (local == false);
		if (fx_candidate < best) {
			best = fx_candidate;
			bits_to_double(candidate, sol, LO, HI, bit_length, total_length);
		}
	}
	f << best << '\n';
	for (auto it : sol)
		f << it << ' ';
	f << '\n';
	auto stop = steady_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	f << "Runtime: " << float(duration.count() / 1000000.00) << "\n\n\n";
}



void Simulated_Annealing(ld(*func)(vector<ld>& x, ll d), ll dim,
	ll ITERATIONS, ld LO, ld HI, const char* file_name)
{

	ld best = numeric_limits<ld>::max();
	ll bit_length = ceil(log2(pow(10, precision) * (HI - LO)));
	ll total_length = bit_length * dim;
	vector<ld> doubles(dim), sol(dim);
	ofstream f(file_name, std::ios_base::app);
	std::random_device rd;
	std::mt19937 gen(rd());
	uniform_real_distribution<> dis(0, 1);

	f << setprecision(10) << fixed;
	constexpr double alpha = 0.7;
	auto start = steady_clock::now();
	for (ll it = 0; it < ITERATIONS; ++it) {
		ld temperatura = 10000;
		cout << it << '\n';
		//int t = 0;
		bitset<BIT_LEN> candidate;
		random_bit_string(candidate, total_length);
		ld fx_candidate = evaluate(candidate, func, dim, LO, HI, bit_length, total_length);
		bool no_downgrade = false;
		if (fx_candidate < best)
			best = fx_candidate;
		do {
			int worse = 50;
			do {
				bitset<BIT_LEN> neighbour = first_improvement(candidate, func, dim, doubles, LO, HI, bit_length, total_length);
				if (neighbour == candidate)
				{
					bitset<BIT_LEN> n2 = candidate;
					ld dif = numeric_limits<ld>::max();
					for (int i = 0; i < total_length; ++i)
					{
						neighbour[i] = !neighbour[i];
						ld fx_neibghour = evaluate(neighbour, func, dim, LO, HI, bit_length, total_length);
						if (fx_candidate - fx_neibghour < dif)
						{
							dif = fx_candidate - fx_neibghour;
							n2 = neighbour;
						}
						neighbour[i] = !neighbour[i];
					}
					neighbour = n2;
				}
				//if (neighbour == candidate)
					//cout << "UAAAAAAAAAAAAAAAAAAAAa\n";
				if (fx_candidate < best)
					best = fx_candidate;

				ld fx_neighbour = evaluate(neighbour, func, dim, LO, HI, bit_length, total_length);

				if (fx_neighbour < fx_candidate) {
					candidate = neighbour;
					fx_candidate = fx_neighbour;
					//cout << "IMPROVE\n";
				}
				else if (dis(gen) < exp(-abs(fx_neighbour - fx_candidate) / temperatura))
				{
					worse--;
					//cout << temperatura << " IES\n";
					if (fx_candidate < best)
						best = fx_candidate;
					candidate = neighbour;
					fx_candidate = fx_neighbour;

				}
				else
					no_downgrade = true;

			} while (worse > 0 && !no_downgrade);
			temperatura = temperatura * alpha;
			//t++;

		} while (temperatura > 0.000000019);
	}
	f << best << '\n';
	for (auto it : sol)
		f << it << ' ';
	f << '\n';
	auto stop = steady_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	f << "Runtime: " << float(duration.count() / 1000000.00) << "\n\n\n";

}
int main() {
	srand(time(NULL));
	//int nr_rulari = 30;

	/// SA 5 dim
	/*
	for (int i = 1; i <= 30; ++i)
	{

		Simulated_Annealing(De_Jongs, 5, 4500, -5.12, 5.12, "SA_de_jong_5.txt");
		Simulated_Annealing(Rastrigin, 5, 4500, -5.12, 5.12, "SA_rastrigin_5.txt");
		Simulated_Annealing(Michalewicz, 5, 4500, 0, pi, "SA_michalewicz_5.txt");
		Simulated_Annealing(Schwefels, 5, 4500, -500, 500, "SA_schwefels_5.txt");
	}
	*/
	/// SA 10 dim
	/*for (int i = 1; i <= 10; ++i)
	{
		Simulated_Annealing(De_Jongs, 10, 1100, -5.12, 5.12, "SA_de_jong_10.txt");
		Simulated_Annealing(Rastrigin, 10, 1100, -5.12, 5.12, "SA_rastrigin_10.txt");
		Simulated_Annealing(Schwefels, 10, 1100, -500, 500, "SA_schwefels_10.txt");
	}*/
	//Simulated_Annealing(Michalewicz, 10, 1100, 0, pi, "SA_michalewicz_10.txt");

	/// SA 30 dim

	/*for (int i = 1; i <= 10; ++i)
	{
		Simulated_Annealing(De_Jongs, 30, 300, -5.12, 5.12, "SA_de_jong_30.txt");
		Simulated_Annealing(Rastrigin, 30, 300, -5.12, 5.12, "SA_rastrigin_30.txt");
		Simulated_Annealing(Michalewicz, 30, 300, 0, pi, "SA_michalewicz_30.txt");
		Simulated_Annealing(Schwefels, 30, 300, -500, 500, "SA_schwefels_30.txt");
	}*/
	/*
	// ---------------- De Jongs -----------------
		Hill_Climber(De_Jongs, worst_improvement, 5, 10000, -5.12, 5.12, "de_jong_worst_5.txt");
	for (int i = 1; i <= nr_rulari; ++i)
		Hill_Climber(De_Jongs, worst_improvement, 10, 5000, -5.12, 5.12, "de_jong_worst_10.txt");
		*/
		//Hill_Climber(De_Jongs, worst_improvement, 30, 500, -5.12, 5.12, "de_jong_worst_30.txt");
// ---------------- Rastrigin -----------------
//run(30, Hill_Climber(Rastrigin, best_improvement, 5, 10000, -5.12, 5.12, "rastrigin_best_5.txt"));
//run(30, Hill_Climber(Rastrigin, first_improvement, 5, 10000, -5.12, 5.12, "rastrigin_first_5.txt"));
//run(30, Hill_Climber(Rastrigin, worst_improvement, 5, 10000, -5.12, 5.12, "rastrigin_worst_5.txt"));

//run(30,Hill_Climber(Rastrigin, best_improvement, 10, 5000, -5.12, 5.12, "rastrigin_best_10.txt"));
//run(30,Hill_Climber(Rastrigin, first_improvement, 10, 5000, -5.12, 5.12, "rastrigin_first_10.txt"));
//run(30,Hill_Climber(Rastrigin, worst_improvement, 10, 5000, -5.12, 5.12, "rastrigin_worst_10.txt"));

//run(30,Hill_Climber(Rastrigin, best_improvement, 30, 500, -5.12, 5.12, "rastrigin_best_30.txt"));
//run(30,Hill_Climber(Rastrigin, first_improvement, 30, 500, -5.12, 5.12, "rastrigin_first_30.txt"));
//run(30, Hill_Climber(Rastrigin, worst_improvement, 30, 500, -5.12, 5.12, "rastrigin_worst_30.txt"));



// ---------------- Michalewicz -----------------
//run(27, Hill_Climber(Michalewicz, best_improvement, 5, 10000, 0, pi, "michalewicz_best_5.txt"));
//run(30, Hill_Climber(Michalewicz, first_improvement, 5, 10000, 0, pi, "michalewicz_first_5.txt"));
//run(30, Hill_Climber(Michalewicz, worst_improvement, 5, 10000, 0, pi, "michalewicz_worst_5.txt"));

//run(30, Hill_Climber(Michalewicz, best_improvement, 10, 5000, 0, pi, "michalewicz_best_10.txt"));
//run(30, Hill_Climber(Michalewicz, first_improvement, 10, 5000, 0, pi, "michalewicz_first_10.txt"));
//run(30, Hill_Climber(Michalewicz, worst_improvement, 10, 5000, 0, pi, "michalewicz_worst_10.txt"));

//run(30, Hill_Climber(Michalewicz, best_improvement, 30, 500, 0, pi, "michalewicz_best_30.txt"));
//run(30, Hill_Climber(Michalewicz, first_improvement, 30, 500, 0, pi, "michalewicz_first_30.txt"));
//run(30, Hill_Climber(Michalewicz, worst_improvement, 30, 500, 0, pi, "michalewicz_worst_30.txt"));


// ---------------- Schwefels -----------------
//run(30, Hill_Climber(Schwefels, best_improvement, 5, 20000, -500, 500, "schwefels_best_5.txt"));
//run(30, Hill_Climber(Schwefels, first_improvement, 5, 20000, -500, 500, "schwefels_first_5.txt"));
//run(30, Hill_Climber(Schwefels, worst_improvement, 5, 20000, -500, 500, "schwefels_worst_5.txt"));

//run(30,Hill_Climber(Schwefels, best_improvement, 10, 10000, -500, 500, "schwefels_best_10.txt"));
//run(30,Hill_Climber(Schwefels, first_improvement, 10, 10000, -500, 500, "schwefels_first_10.txt"));
//run(30,Hill_Climber(Schwefels, worst_improvement, 10, 10000, -500, 500, "schwefels_worst_10.txt"));

//run(30,Hill_Climber(Schwefels, best_improvement, 30, 500, -500, 500, "schwefels_best_30.txt"));
//run(30,Hill_Climber(Schwefels, first_improvement, 30, 500, -500, 500, "schwefels_first_30.txt"));
//run(14,Hill_Climber(Schwefels, worst_improvement, 30, 500, -500, 500, "schwefels_worst_30.txt")); 
	return 0;
}
