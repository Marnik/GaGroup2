package cifo;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Random;

public class GeneticAlgorithm extends SearchMethod {

	protected ProblemInstance instance;
	protected int populationSize, numberOfGenerations;
	protected double mutationProbability;
	protected int tournamentSize;
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
		tournamentSize = Main.TOURNAMENT_SIZE;
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
			Solution[] offspring = new Solution[populationSize];
			for (int k = 0; k < population.length; k++) {
				//int[] parents = selectParents();
				int[] parents = rouletteSelection();

				offspring[k] = applyCrossover(parents);
				if (r.nextDouble() <= mutationProbability) {
					offspring[k] = offspring[k].applyMutation();
				}
				offspring[k].evaluate();
			}

			population = survivorSelection(offspring);
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
			// the individual containing the range in which the random number is located.
			for (int j = 0; j < populationSize; j++) { 
				if (randomDouble < rouletteWheel.get(j)) {
					parents[i] = j;
					break;
				}
			}
		}
		return parents;
	}
	
	// Method that is used to set up the roulette wheel which is used to select parents
	private HashMap<Integer, Double> createRouletteWheel(){
		double totalFitness = 0.0;
		HashMap<Integer, Double> invertedProbabilities = new HashMap<Integer, Double>();
		HashMap<Integer, Double> rouletteWheel = new HashMap<Integer, Double>();
		
		for (int i = 0; i < populationSize; i++) {
			totalFitness += population[i].getFitness();
		}

		// Minimization problem, so invert values that are normally used for a maximization problem.
		for (int i = 0; i < populationSize; i++) {
			invertedProbabilities.put(i, 1 - (population[i].getFitness() / totalFitness));
		}

		// We want to repeat the process for found values without
		// subtracting from 1. This time we'll assign the cumulative probabilities.
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

	public Solution applyCrossover(int[] parents) {
		Solution firstParent = population[parents[0]];
		Solution secondParent = population[parents[1]];
		Solution offspring = firstParent.copy();
		int crossoverPoint = r.nextInt(instance.getNumberOfTriangles() * Solution.VALUES_PER_TRIANGLE);
		for (int i = crossoverPoint; i < instance.getNumberOfTriangles() * Solution.VALUES_PER_TRIANGLE; i++) {
			offspring.setValue(i, secondParent.getValue(i));
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
	
	private double getTotalFitness(){
		double totalFitness = 0.0;
		HashMap<Integer, Double> invertedProbabilities = new HashMap<Integer, Double>();
		
		for (int i = 0; i < populationSize; i++) {
			totalFitness += population[i].getFitness();
		}

		// Minimization problem, so invert values that are normally used for a maximization problem.
		for (int i = 0; i < populationSize; i++) {
			invertedProbabilities.put(i, 1 - (population[i].getFitness() / totalFitness));
		}

		// We want to repeat the process for found values without
		// subtracting from 1. This time we'll assign the cumulative probabilities.
		totalFitness = 0.0;
		for (double d : invertedProbabilities.values()) {
			totalFitness += d;
		}
		
		return totalFitness;
	}
}
