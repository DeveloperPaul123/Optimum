#ifndef GENETICALGORITHM_H
#define GENETICALGORITHM_H
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>

namespace Optimum {
	/**
	* Class that encapsulates operations with a Chromosome in the Genetic Algorithm.
	*/
	template <typename T> class Chromosome{
		std::vector<T> mValues;
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
		Chromosome(T* values, int size, int mutationSize = 2) {
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
		* [] Operator to make looking up values easy.
		* @param i the index to look up.
		* @return T the value at the given index.
		*/
		T& operator[](int i) {
			return mValues[i];
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
			int position = getRandomInt(0, mValues.size() - 1);
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
			Chromosome<T> chrome;
			chrome.prepare(other.size());
			for (int i = 0; i < other.size(); i++) {
				T value = this->at(i) + other.at(i);
				T avg = value / 2;
				chrome.set(avg, i);
			}

			return chrome;
		}

		/**
		* Performs an addition cross over.
		* @param other the other chromosome to cross over with.
		*/
		Chromosome<T> addCrossOver(Chromosome<T> other) {
			Chromosome<T> chrome;
			chrome.prepare(other.size());
			for (int i = 0; i < other.size(); i++) {
				chrome.set(this->at(i) + other.at(i), i);
			}

			return chrome;
		}

		/**
		* Subtractiong cross over. Subtracts values from the two chromosomes.
		* @param other the other chromosome.
		*/
		Chromosome<T> subtractCrossOver(Chromosome<T> other) {
			Chromosome<T> chrome;
			chrome.prepare(other.size());
			for (int i = 0; i < other.size(); i++) {
				chrome.set(this->at(i) - other.at(i), i);
			}
			return chrome;
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
			srand(static_cast <unsigned> (time(0)));
			for (int i = 0; i < mutationSize; i++) {
				int pos = getRandomInt(0, size());
				if (pos > size() - 1) {
					pos -= 1;
				}
				T value = at(pos);
				T newValue = (value*-1);
				chrome.set(newValue, pos);
			}
			return chrome;
		}

		/**
		* Creates a random mutation on this Chromosome and returns
		* a newly mutated chromosome.
		* @return Chromosome<T> the newly mutated chromosome.
		*/
		Chromosome<T> randomMutation() {
			Chromosome<T> chrome = *this;
			srand(static_cast <unsigned> (time(0)));
			for (int i = 0; i < size(); i++) {
				T value = at(i);
				T newValue = value * getRandomFloat(0.1f, 2.0f);
				chrome.set(newValue, i);
			}
			return chrome;
		}

		/**
		* Helper method that returns a random float between a given range.
		* @param LO the minimum number.
		* @param HI the maximum number.
		*/
		float getRandomFloat(float LO, float HI) {
			float r3 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
			return r3;

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

	/**
	* This class encapsualtes a population of chromosomes from which to mutate, select and do
	* other operations to.
	*/
	template <class T> class Population {

		std::vector<Chromosome<T>> chromes;
		int popSize;
		float eliteRate, mutRate;

	public:

		/**
		* Default constructor for a population. Holds many members (chromosomes).
		* @param popSize the size of the population, default is 1000.
		* @param eliteRate the elitist rate of the population, default is 10%.
		* @param mutationRate the mutation rate of the population, default is 40%.
		*/
		Population(int popSize = 1000, float eliteRate = 0.10, float mutationRate = 0.40) {
			srand(static_cast <unsigned> (time(0)));
			this->popSize = popSize;
			this->eliteRate = eliteRate;
			this->mutRate = mutationRate;
		}

		/**
		* Performs a bubble sort on the population. Sorts it in ascending order.
		*/
		void bubbleSort() {
			int i, j, flag = 1;    // set flag to 1 to start first pass
			Chromosome<T> temp;             // holding variable
			int numLength = chromes.size();
			for (i = 1; (i <= numLength) && flag; i++)
			{
				flag = 0;
				for (j = 0; j < (numLength - 1); j++)
				{
					if (chromes[j + 1].getFitness() < chromes[j].getFitness()) {
						// swap elements
						temp = chromes[j];
						temp.setFitness(chromes[j].getFitness());
						chromes[j] = chromes[j + 1];
						chromes[j + 1] = temp;
						// indicates that a swap occurred.
						flag = 1;
					}
				}
			}
		}

		/**
		* Adds a member to this population.
		* @param member the new member to add.
		*/
		void addMember(Chromosome<T> member) {
			chromes.push_back(member);
		}

		/**
		* Returns a member of this population at a given position.
		* @param pos the position to look up.
		*/
		Chromosome<T> getMember(int pos) {
			return chromes[pos];
		}

		/**
		* Set the values of this population equal to the members of another given population.
		* Basically copies all the values from one population into this one.
		* @param other the population to copy.
		*/
		void set(Population<T> other) {
			chromes.clear();
			for (int i = 0; i < other.size(); i++) {
				chromes.push_back(other.getMember(i));
			}
		}

		/**
		* Replaces a member of this population with a new member.
		* @param pos the position of the new member.
		* @param member the new member.
		*/
		void setMember(int pos, Chromosome<T> member) {
			chromes[pos] = member;
		}

		/**
		* Performs elitism selection. Returns a population with the best members of this
		* current population. Nothing is done to these "best" candidates.
		* @return Population<T> the "elite" population.
		*/
		Population<T> elitism() {
			int ne = eliteRate * popSize;
			bubbleSort();
			Population<T> pop(this->popSize, this->eliteRate, this->mutRate);
			for (int i = 0; i < ne; i++) {
				pop.addMember(getMember(i));
			}

			return pop;
		}

		/**
		* Generates a crossover population from the current population. Does this by
		* randomly choosing two members of the current population as parents and then
		* mixing their values across two children. Puts the children in a new population.
		* The parents are not altered in any way.
		* @return Population<T> the new population with the bred children.
		*/
		Population<T> generateCrossoverPopulation() {
			int ne = eliteRate*popSize;
			int nc = (popSize - ne) / 2;
			bubbleSort();
			Population<T> newPop(this->popSize, this->eliteRate, this->mutRate);

			//TODO: Roulette Wheel Selection or Tournament Selection. 

			for (int i = 0; i < nc; i++) {

				int pos1 = getRandomInt(0, chromes.size());
				int pos2 = getRandomInt(0, chromes.size());
				if (pos1 == pos2) {
					pos2 = getRandomInt(0, chromes.size());
				}

				int pos3 = getRandomInt(0, chromes.size());
				int pos4 = getRandomInt(0, chromes.size());
				if (pos4 == pos3) {
					pos4 = getRandomInt(0, chromes.size());
				}

				Chromosome<T> first = chromes[pos1];
				Chromosome<T> second = chromes[pos2];
				Chromosome<T> third = chromes[pos3];
				Chromosome<T> fourth = chromes[pos4];
				Chromosome<T> crossOne;
				Chromosome<T> crossTwo;
				Chromosome<T> crossThree;
				Chromosome<T> crossFour;

				crossOne.prepare(first.size());
				crossTwo.prepare(second.size());
				crossThree.prepare(third.size());
				crossFour.prepare(fourth.size());

				for (int j = 0; j < first.size(); j++) {
					float beta = getRandomFloat(0.0f, 1.0f);
					T value = beta* first.at(j) + (1 - beta)*second.at(j);
					T valueTwo = beta*third.at(j) + (1 - beta)*fourth.at(j);
					crossOne.set(value, j);
					crossThree.set(valueTwo, j);
				}


				for (int n = 0; n < first.size(); n++) {
					float beta = getRandomFloat(0.0f, 1.0f);
					T value = beta* first.at(n) + (1 - beta)*second.at(n);
					T valueTwo = beta*third.at(n) + (1 - beta)*fourth.at(n);
					crossTwo.set(value, n);
					crossFour.set(value, n);
				}

				newPop.addMember(crossOne);
				newPop.addMember(crossTwo);
			}

			return newPop;
		}

		/**
		* Adds a given population to this population.
		* @param pop the population to add.
		* @return Population<T> the new population with the members from this
		*		population and the other population (i.e. pop).
		*/
		Population<T> addPopulation(Population<T> pop) {
			Population<T> newPop;
			for (int i = 0; i < size(); i++) {
				newPop.addMember(getMember(i));
			}
			for (int i = 0; i < pop.size(); i++) {
				newPop.addMember(pop.getMember(i));
			}

			return newPop;
		}

		/**
		* Returns the size of this population.
		* @return int the size.
		*/
		int size() {
			return chromes.size();
		}

	private:

		/**
		* Helper method that returns a random float between a given range.
		* @param LO the minimum.
		* @param HI the maximum.
		* @return a random number between LO and HI.
		*/
		float getRandomFloat(float LO, float HI) {
			float r3 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
			return r3;

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

	/**
	* Genetic Algorithm Class. Holds the logic for stepping through generations and
	* handling the new populations.
	*/
	template <class T> class GeneticAlgorithm {
		int popSize;
		int maxIter;
		float elite;
		float mutation;
		int iter;
		Population<T> pop;

	public:

		/**
		* Constructor. Initializes various values.
		* @param popSize the population size, the default is 200.
		* @param maxIter the maximum number of generations to create. The default is 10,000.
		* @param elite the elitest rate, default it 10%.
		* @param mutation rate, the default is 25%.
		*/
		GeneticAlgorithm(int popSize = 200, int maxIter = 10000, float elite = 0.10f, float mutation = 0.25f) {
			srand(static_cast <unsigned> (time(0)));
			this->popSize = popSize;
			this->maxIter = maxIter;
			this->elite = elite;
			this->mutation = mutation;
			this->iter = 0;
			pop = Population<T>(popSize, elite, mutation);
		}

		/**
		* Checks to see if the algorithm is done. Currently only looks if
		* the number of iterations has passed or is equal to the max number
		* of iterations.
		* @return bool true if done, false otherwise.
		*/
		bool done() {
			if (iter >= maxIter) {
				return true;
			}
			else {
				return false;
			}
		}

		/**
		* Generates a random population given an "initial guess".
		* @param initialGuess the initial guess.
		* @return Population<float> a population of size popSize.
		*/
		Population<float> generatePopulation(Chromosome<float> initialGuess) {
			for (int i = 0; i < popSize; i++) {
				Chromosome<float> chrome;
				chrome.prepare(initialGuess.size());
				for (int i = 0; i < initialGuess.size(); i++) {
					float val = getRandomFloat(initialGuess.at(i) - 20.0f, initialGuess.at(i) + 20.0f);
					chrome.set(val, i);
				}
				pop.addMember(chrome);
			}
			return pop;
		}

		/**
		* Generates an random population given an initial integer guess.
		* @param initial guess.
		* @return Population<int> a newly generated population of size popSize.
		*/
		Population<int> generatePopulation(Chromosome<int> initialGuess) {

			for (int i = 0; i < popSize; i++) {
				Chromosome<int> chrome;
				chrome.prepare(initialGuess.size());
				for (int i = 0; i < initialGuess.size(); i++) {
					int val = getRandomInt(initialGuess.at(i) - 500, initialGuess.at(i) + 500);
					chrome.set(val, i);
				}
				pop.addMember(chrome);
			}

			return pop;
		}

		/**
		* Take a step with the population and generate a new generation.
		* Goes through the steps of elitism, crossing over, mutation and updating the current
		* population.
		* @param pop the current population.
		* @return the new population.
		*/
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
			pop.set(pop1.addPopulation(pop2));
			pop.bubbleSort();

			iter += 1;

			return pop;
		}

		/**
		* Helper method that generates a random float in a given range.
		* @param LO the minimum number.
		* @param HI the maximum number.
		*/
		float getRandomFloat(float LO, float HI) {
			float r3 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
			return r3;
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
}
#endif GENETICALGORITHM_H