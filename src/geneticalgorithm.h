#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>


using namespace std;

/**
* Class that encapsulates operations with a Chromosome in the Genetic Algorithm. 
*/
template <typename T> class Chromosome{
	vector<T> mValues;
protected:
	float fitness;
	int mutationSize;
public:

	/**
	* Default constructor, sets the mutation size to 2 as a default. 
	* @param mutationSize the max number of mutations and the max size of the mutation. 
	*/
	Chromosome(int mutationSize = 2){
		srand(static_cast <unsigned> (time(0)));
		this->mutationSize = mutationSize;
		this->fitness = 0.0;
	}

	/**
	* Constructor with initialization of values, the size of the values and the mutation size.
	* @param values the array of values for this Chromosome. 
	* @param size the size of the array of values. 
	* @param mutationSize the max size of mutations, default is 2. 
	*/
	Chromosome(T* values, int size, int mutationSize=2) {
		srand(static_cast <unsigned> (time(0)));
		for (int i = 0; i < size; i++) {
			mValues.push_back(values[i]);
		}
		this->mutationSize = mutationSize;
		this->fitness = 0.0;
	}

	/**
	* Prepare the data in this chromosome. 
	* @param size the size of the incoming data.
	*/
	void prepare(int size) {
		for (int i = 0; i < size; i++) {
			mValues.push_back(0);
		}
	}

	/**
	* Set the value of a given position in the chromosome. 
	* @param value class T value. 
	* @param pos the position to place the new value at. 
	*/
	void set(T value, int pos) {
		mValues[pos] = value;
	}

	/**
	* Read the value of the chromosome at a given point. 
	* @param pos the position to read from.
	* @return T the value. 
	*/
	T at(int pos) {
		return mValues[pos];
	}

	/**
	* Gets the size of this chromosome. 
	* @return int the size of data in the chromosome. 
	*/
	int size() {
		return mValues.size();
	}

	/**
	* Set the fitness of this chromosome. 
	* @param fitness the fitness of this chromosome. 
	*/
	void setFitness(float fitness) { this->fitness = fitness; }

	/**
	* Read the fitness of this chromosome. 
	* @return float the fitness of this chromosome. 
	*/
	float getFitness() { return this->fitness; }
 
	/**
	* Performs a single point crossover of the chromosome at a random index. 
	* @param other the other Chromosome to cross over with. 
	* @return Chromosome<T> a new Chromosome. 
	*/
	Chromosome<T> singlePointCrossOver(Chromosome<T> other) {
		int position = getRandomInt(0, mValues.size()-1);
		Chromosome<T> chrome;
		chrome.prepare(other.size());
		for (int i = 0; i < position; i++) {
			chrome.set(at(i), i);
		}
		for (int j = position; j < other.size(); j++) {
			chrome.set(other.at(j), j);
		}
		return chrome;
	}

	/**
	* Performs a double point cross over.
	* @param other the chromosome to cross over with. 
	* @return Chromosome<T> the new Chromosome<T>.
	*/
	Chromosome<T> doublePointCrossOver(Chromosome<T> other) {
		int posOne = getRandomInt(0, mValues.size() - 1);
		int posTwo = getRandomInt(0, mValues.size() - 1);
		while (posTwo <= posOne) {
			if (posOne == mValues.size() - 1) {
				posOne = getRandomInt(0, mValues.size() - 1);
			}
			posTwo = getRandomInt(0, mValues.size() - 1);
		}

		Chromosome<T> chrome;
		chrome.prepare(other.size());
		for (int i = 0; i < posOne; i++) {
			chrome.set(at(i), i);
		}

		for (int j = posOne; j < posTwo; j++) {
			chrome.set(other.at(j), j);
		}

		for (int n = posTwo; n < mValues.size(); n++) {
			chrome.set(at(n), n);
		}
		return chrome;
	}

	/**
	* Performs a midpoint cross over. Basically averages the two values. 
	* @param other the other chromosome to cross over with. 
	*/
	Chromosome<T> midPointCrossOver(Chromosome<T> other) {

	}

	/**
	* Performs an addition cross over. 
	* @param other the other chromosome to cross over with. 
	*/
	Chromosome<T> addCrossOver(Chromosome<T> other) {

	}

	/**
	* Subtractiong cross over. Subtracts values from the two chromosomes. 
	* @param other the other chromosome. 
	*/
	Chromosome<T> subtractCrossOver(Chromosome<T> other) {

	}
	
	/**
	* Swap mutation, swaps multiple positions in the chromosome with new values. 
	* @return Chromosome<T> the new, mutated chromosome. 
	*/
	Chromosome<T> swapMutation(){
		Chromosome<T> chrome = *this;
	}

	/**
	* Performs a inversion mutation on this chromosome. Inverts the number
	* by multiplying it by negative one. 
	* @return Chromosome<T> the new chromosome. 
	*/
	Chromosome<T> invertMutation() {
		Chromosome<T> chrome = *this;
		for (int i = 0; i < mutationSize; i++) {
			int pos = getRandomInt(0, size() - 1);
			T value = at(i);
			T newValue = value*-1;
			chrome.set(pow, newValue);
		}
		return chrome;
	}

	/**
	* Helper method that returns a random integer between a given range. 
	* @param min the minimum number.
	* @param max the maximum number. 
	*/
	int getRandomInt(int min, int max) {
		return rand() % max + min;
	}
};


template <class T> class Population {

	vector<Chromosome<T>> chromes;
	int popSize, maxIter;
	float eliteRate, mutRate;

public:
	Population(int popSize = 1000, int maxIterations = 100000, float eliteRate = 0.10, float mutationRate = 0.25f) {
		srand(static_cast <unsigned> (time(0)));
		this->popSize = popSize;
		this->maxIter = maxIterations;
		this->eliteRate = eliteRate;
		this->mutRate = mutationRate;
	}

	void addMember(Chromosome<T> member) {
		chromes.push_back(member);
	}
	Chromosome<T> getMember(int pos) {
		return chromes[pos];
	}
	void setMember(int pos, Chromosome<T> member) {
		chroms[pos] = member;
	}
	
	Population<T> elitism() {
		int ne = eliteRate * popSize;
		sortPopulation();
		Population<T> pop(this->popSize, this->maxIter, this->eliteRate, this->mutRate);
		for (int i = 0; i < ne; i++) {
			pop.addMember(getMember(i));
		}

		return pop;
	}

	Population<T> generateCrossoverPopulation() {
		int ne = eliteRate*popSize;
		int nc = (popSize - ne) / 2;

		Population<T> newPop(this->popSize, this->maxIter, this->eliteRate, this->mutRate);

		for (int i = 0; i < nc) {
			int pos1 = getRandomInt(0, chromes.size());
			int pos2 = getRandomInt(0, chromes.size());
			while (pos1 == pos2) {
				int pos1 = getRandomInt(0, chromes.size());
				int pos2 = getRandomInt(0, chromes.size());
			}

			Chromosome<T> first = chromes[pos1];
			Chromosome<T> second = chromes[pos2];
			Chromosome<T> crossOne = first.singlePointCrossOver(second);
			Chromosome<T> crossTwo = second.singlePointCrossOver(first);
			newPop.addMember(crossOne);
			newPop.addMember(crossTwo);
		}
	}

	void addPopulation(Population<T> pop) {
		for (int i = 0; i < pop.size(); i++) {
			this->chromes.push_back(pop.getMember(i));
		}
	}

	void sortPopulation() {
		sort(chromes.begin(), chromes.end(), chromSort);
	}

	int size() {
		return popSize;
	}
private:
	bool chromSort(Chromosome<T> chrom1, Chromosome<T> chrom2) {
		return (chrom1.getFitness() < chrom2.getFitness());
	}

	/**
	* Helper method that returns a random integer between a given range.
	* @param min the minimum number.
	* @param max the maximum number.
	*/
	int getRandomInt(int min, int max) {
		return rand() % max + min;
	}
};

template <class T> class GeneticAlgorithm {
	int popSize;
	int maxIter;
	float elite;
	float mutation;
	Population<T> pop;
public:
	GeneticAlgorithm(int popSize = 1000, int maxIter = 10000, float elite = 0.10f, float mutation = 0.25f) {
		srand(static_cast <unsigned> (time(0)));
		this->popSize = popSize;
		this->maxIter = maxIter;
		this->elite = elite;
		this->mutation = mutation;
	}

	void generatePopulation(Chromosome<float> initialGuess) {
		for (int i = 0; i < popSize; i++) {
			Chromosome<float> chrome;
			chrome.prepare(initialGuess.size());
			for (int i = 0; i < initialGuess.size(); i++) {
				float val = getRandomFloat(initialGuess.at(i) - 500.0f, initialGuess.at(i) + 500.0f);
				chrome.set(val, i);
			}
			pop.addMember(chrome);
		}
	}

	void generatePopulation(Chromosome<int> initialGuess) {
	
		for (int i = 0; i < popSize; i++) {
			Chromosome<int> chrome;
			chrome.prepare(initialGuess.size());
			for (int i = 0; i < initialGuess.size(); i++) {
				int val = getRandomInt(initialGuess.at(i) - 500, initialGuess.at(i) + 500);
				chrome.set(val, i);
			}
			pop.addMember(chrome);
		}
	}

	Population<T> step(Population<T> pop) {
		//elitism. 
		Population<T> pop1 = pop.elitism();
		//cross over. 
		Population<T> pop2 = pop.generateCrossoverPopulation();

		int ne = elite * popSize;
		int nc = (popSize - ne) / 2;

		//mutation. 
		for (int i = 0; i < nc; i++) {
			Chromosome<T> chrome = pop2.getMember(i);
			Chromosome<T> newChrome = chrome.invertMutation();
			pop2.setMember(i, newChrome);
		}

		//update.
		this->pop = pop1.addPopulation(pop2);
		return this->pop;
	}

	float getRandomFloat(float LO, float HI) {
		float r3 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
		return r3. 
	}

	/**
	* Helper method that returns a random integer between a given range.
	* @param min the minimum number.
	* @param max the maximum number.
	*/
	int getRandomInt(int min, int max) {
		return rand() % max + min;
	}
};


#endif GENETICALGORITHM_H