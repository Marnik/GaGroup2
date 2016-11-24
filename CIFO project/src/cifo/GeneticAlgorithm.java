package cifo;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Random;

public class GeneticAlgorithm extends SearchMethod {

	protected ProblemInstance instance;
	protected int populationSize, numberOfGenerations;
	protected double mutationProbability;
	protected int tournamentSize;
	protected int numberOfOffsprings;
	protected int numberOfTriangleMutations;
	protected boolean printFlag;
	protected Solution currentBest;
	protected int currentGeneration;
	protected Solution[] population;
	protected Random r;

	public GeneticAlgorithm() {
		instance = new ProblemInstance(Main.NUMBER_OF_TRIANGLES);
		populationSize = Main.POPULATION_SIZE;
		numberOfGenerations = Main.NUMBER_OF_GENERATIONS;
		mutationProbability = Main.MUTATION_PROBABILIY;
		numberOfTriangleMutations = Main.NUMBER_OF_TRIANGLE_MUTATIONS;
		tournamentSize = Main.TOURNAMENT_SIZE;
		numberOfOffsprings = Main.NUMBER_OF_OFFSPRINGS;
		printFlag = false;
		currentGeneration = 0;
		r = new Random();
	}

	public void run() {
		double startTime = System.currentTimeMillis();
		initialize();
		evolve();
		double endTime = System.currentTimeMillis();
		double executionTime = (endTime - startTime) / 60000;
		Main.addBestSolution(currentBest, executionTime);
	}

	public void initialize() {
		population = new Solution[populationSize];
		for (int i = 0; i < population.length; i++) {
			population[i] = new Solution(instance);
			population[i].evaluate();
		}
		updateCurrentBest();
		updateInfo();
		currentGeneration++;
	}

	public void updateCurrentBest() {
		currentBest = getBest(population);
	}

	public void evolve() {
		while (currentGeneration <= numberOfGenerations) {
			Solution[] offspring = new Solution[numberOfOffsprings];
			Solution[] nextPopulation = new Solution[populationSize];
			int nextPopulationSize = 0;

			// Used to choose random triangles to mutate
			int[] triangleIndexes = new int[instance.getNumberOfTriangles()];
			for (int t = 0; t < instance.getNumberOfTriangles(); t++) {
				triangleIndexes[t] = t;
			}

			for (int k = 0; nextPopulationSize < population.length; k += numberOfOffsprings) {
				int[] parents = selectParents();
				// int[] parents = rouletteSelection();
				
				offspring = applyCrossover(parents);
				triangleIndexes = shuffleArray(triangleIndexes);
				
				// Loop through either one or two offsprings
				for (int j = 0; j < numberOfOffsprings; j++) {
					if (r.nextDouble() <= mutationProbability) {
						for (int tm = 0; tm < numberOfTriangleMutations; tm++) {
							offspring[j] = offspring[j].applyMutation(triangleIndexes[tm]);
						}
					}
					offspring[j].evaluate();

					// Generating two offsprings resulted in one extra offspring
					// if this converts to true.
					// Eliminate the offspring with the lowest fitness
					if (nextPopulationSize == population.length) {
						int worstFitnessIndex = getWorstIndex(nextPopulation);
						if (offspring[j].getFitness() < nextPopulation[worstFitnessIndex].getFitness()) {
							nextPopulation[worstFitnessIndex] = offspring[j];

						}
					} else {
						nextPopulation[k + j] = offspring[j];
						nextPopulationSize++;
					}

				}
			}

			population = survivorSelection(nextPopulation);
			updateCurrentBest();
			updateInfo();
			currentGeneration++;
		}
	}

	// Method to select parents using the tournament selection method
	public int[] selectParents() {
		int[] parents = new int[2];
		parents[0] = r.nextInt(populationSize);
		for (int i = 0; i < tournamentSize; i++) {
			int temp = r.nextInt(populationSize);
			if (population[temp].getFitness() < population[parents[0]].getFitness()) {
				parents[0] = temp;
			}
		}

		parents[1] = r.nextInt(populationSize);
		for (int i = 0; i < tournamentSize; i++) {
			int temp = r.nextInt(populationSize);
			if (population[temp].getFitness() < population[parents[1]].getFitness()) {
				parents[1] = temp;
			}
		}
		return parents;
	}

	// Method to select parents using the roulette wheel selection method
	private int[] rouletteSelection() {
		int[] parents = new int[2];
		HashMap<Integer, Double> rouletteWheel = createRouletteWheel();

		// 'Spin the roulette wheel twice' so select two parents
		for (int i = 0; i < 2; i++) {
			Random r = new Random();
			double randomDouble = 1.0 * r.nextDouble();

			// Loop through the HashMap containing the probabilities and select
			// the individual containing the range in which the random number is
			// located.
			for (int j = 0; j < populationSize; j++) {
				if (randomDouble < rouletteWheel.get(j)) {
					parents[i] = j;
					break;
				}
			}
		}
		return parents;
	}

	// Method that is used to set up the roulette wheel which is used to select
	// parents
	private HashMap<Integer, Double> createRouletteWheel() {
		double totalFitness = 0.0;
		HashMap<Integer, Double> invertedProbabilities = new HashMap<Integer, Double>();
		HashMap<Integer, Double> rouletteWheel = new HashMap<Integer, Double>();

		for (int i = 0; i < populationSize; i++) {
			totalFitness += population[i].getFitness();
		}

		// Minimization problem, so invert values that are normally used for a
		// maximization problem.
		for (int i = 0; i < populationSize; i++) {
			invertedProbabilities.put(i, 1 - (population[i].getFitness() / totalFitness));
		}

		// We want to repeat the process for found values without
		// subtracting from 1. This time we'll assign the cumulative
		// probabilities.
		totalFitness = 0.0;
		for (double d : invertedProbabilities.values()) {
			totalFitness += d;
		}

		for (int i = 0; i < populationSize; i++) {
			if (i == 0) {
				rouletteWheel.put(i, invertedProbabilities.get(i) / totalFitness);
			} else {
				rouletteWheel.put(i, (invertedProbabilities.get(i) / totalFitness) + rouletteWheel.get(i - 1));
			}
		}

		return rouletteWheel;
	}

	public Solution[] applyCrossover(int[] parents) {
		Solution firstParent = population[parents[0]];
		Solution secondParent = population[parents[1]];
		Solution offspring[] = new Solution[numberOfOffsprings];
		int crossoverPoint = r.nextInt(instance.getNumberOfTriangles() * Solution.VALUES_PER_TRIANGLE);
		offspring[0] = firstParent.copy();

		for (int i = crossoverPoint; i < instance.getNumberOfTriangles() * Solution.VALUES_PER_TRIANGLE; i++) {
			offspring[0].setValue(i, secondParent.getValue(i));
		}

		if (numberOfOffsprings == 2) {
			offspring[1] = secondParent.copy();
			for (int i = crossoverPoint; i < instance.getNumberOfTriangles() * Solution.VALUES_PER_TRIANGLE; i++) {
				offspring[1].setValue(i, firstParent.getValue(i));
			}
		}

		return offspring;
	}

	public Solution[] survivorSelection(Solution[] offspring) {
		Solution bestParent = getBest(population);
		Solution bestOffspring = getBest(offspring);
		if (bestOffspring.getFitness() <= bestParent.getFitness()) {
			return offspring;
		} else {
			Solution[] newPopulation = new Solution[population.length];
			newPopulation[0] = bestParent;
			int worstOffspringIndex = getWorstIndex(offspring);
			for (int i = 0; i < newPopulation.length; i++) {
				if (i < worstOffspringIndex) {
					newPopulation[i + 1] = offspring[i];
				} else if (i > worstOffspringIndex) {
					newPopulation[i] = offspring[i];
				}
			}
			return newPopulation;
		}
	}

	public Solution getBest(Solution[] solutions) {
		Solution best = solutions[0];
		for (int i = 1; i < solutions.length; i++) {
			if (solutions[i].getFitness() < best.getFitness()) {
				best = solutions[i];
			}
		}
		return best;
	}

	public int getWorstIndex(Solution[] solutions) {
		Solution worst = solutions[0];
		int index = 0;
		for (int i = 1; i < solutions.length; i++) {
			if (solutions[i].getFitness() > worst.getFitness()) {
				worst = solutions[i];
				index = i;
			}
		}
		return index;
	}

	public void updateInfo() {
		currentBest.draw();
		series.add(currentGeneration, currentBest.getFitness());
		if (printFlag) {
			System.out.printf("Generation: %d\tFitness: %.1f\n", currentGeneration, currentBest.getFitness());
		}
	}

	private double getTotalFitness() {
		double totalFitness = 0.0;
		HashMap<Integer, Double> invertedProbabilities = new HashMap<Integer, Double>();

		for (int i = 0; i < populationSize; i++) {
			totalFitness += population[i].getFitness();
		}

		// Minimization problem, so invert values that are normally used for a
		// maximization problem.
		for (int i = 0; i < populationSize; i++) {
			invertedProbabilities.put(i, 1 - (population[i].getFitness() / totalFitness));
		}

		// We want to repeat the process for found values without
		// subtracting from 1. This time we'll assign the cumulative
		// probabilities.
		totalFitness = 0.0;
		for (double d : invertedProbabilities.values()) {
			totalFitness += d;
		}

		return totalFitness;
	}

	private int[] shuffleArray(int[] ar) {
		// If running on Java 6 or older, use `new Random()` on RHS here
		Random rnd = new Random();
		for (int i = ar.length - 1; i > 0; i--) {
			int index = rnd.nextInt(i + 1);
			// Simple swap
			int a = ar[index];
			ar[index] = ar[i];
			ar[i] = a;
		}

		return ar;
	}
}
