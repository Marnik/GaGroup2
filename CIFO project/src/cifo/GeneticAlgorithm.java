package cifo;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
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
	protected String crossoverMethod;
	protected final String [] mutation_methods = {"standard", "cycle", "PMXO", "PMXO_triangle", "cycle_triangle", "six_way", "seperate", "random_triangle"};
	protected int optimalCount = 20;
	protected int [][] optimalMethodCount;
	protected String [] optimalMethodSelection;

	protected int pType, hammingDistance, maxPossibleDistance;
	protected double stdDev, pVariance, hammingDistanceRatio;

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
		this.crossoverMethod = Main.CROSSOVER_METHOD;
		r = new Random();
		optimalMethodCount = new int [numberOfGenerations/optimalCount][mutation_methods.length];
		optimalMethodSelection = new String [numberOfGenerations/optimalCount];
	}

	public void run() {
		double startTime = System.currentTimeMillis();
		initialize();
		evolve();
		double endTime = System.currentTimeMillis();
		double executionTime = (endTime - startTime) / 60000;
		for (int i=0; i<optimalMethodSelection.length; i++) {
			int mostCommonMethodUsed = 0;
			for (int j=1; j<optimalMethodCount[i].length; j++) {
				if (optimalMethodCount[i][j] > optimalMethodCount[i][mostCommonMethodUsed]) mostCommonMethodUsed = j;
			}
			optimalMethodSelection[i] = mutation_methods[mostCommonMethodUsed];
		}
		Main.addBestSolution(currentBest, executionTime, optimalMethodSelection);
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

			switch (Main.MUTATION_FUNCTION) {
			case "linear":
				numberOfTriangleMutations = getAdaptedLinearTriangleMutation();
				break;
			case "exponential":
				numberOfTriangleMutations = getAdaptedExponentialTriangleMutation();
				break;
			case "prefixed":
				numberOfTriangleMutations = getPrefixedTriangleMutation();
				break;
			}
			
			if (Main.INCREMENT_RATE != 0 && currentGeneration % Main.INCREMENT_RATE == 0) {
				numberOfTriangleMutations++;
			}
			if (Main.DECREMENT_RATE != 0 && numberOfTriangleMutations > 0 && currentGeneration % Main.DECREMENT_RATE == 0) {
				numberOfTriangleMutations--;
			}
			
			Solution[] offspring = new Solution[numberOfOffsprings];
			Solution[] nextPopulation = new Solution[populationSize];
			int nextPopulationSize = 0;

			// Used to choose random triangles to mutate
			int[] triangleIndexes = new int[instance.getNumberOfTriangles()];
			for (int t = 0; t < instance.getNumberOfTriangles(); t++) {
				triangleIndexes[t] = t;
			}

			for (int k = 0; nextPopulationSize < population.length; k += numberOfOffsprings) {

				int[] parents = new int[2];
				if (Main.SELECTION_METHOD == "tournament") {
					parents = tournamentSelection();
				}
				if (Main.SELECTION_METHOD == "roulette") {
					parents = rouletteSelection();
				}

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
					
					if (nextPopulationSize >= population.length) {

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

			generateDiversityMeasures();

			currentGeneration++;
		}
	}

	private int getPrefixedTriangleMutation() {
		int numberOfTriangleMutations = this.numberOfTriangleMutations;
		if (currentGeneration == 200) {
			System.out.println("Changing triangle mutations from " + this.numberOfTriangleMutations + " to 2");
			numberOfTriangleMutations = 2;
		}
		if (currentGeneration == 400) {
			System.out.println("Changing triangle mutations from " + this.numberOfTriangleMutations + " to 1");
			numberOfTriangleMutations = 1;
		}
		return numberOfTriangleMutations;
	}

	private int getAdaptedExponentialTriangleMutation() {
		// f(x) = 8,9885e-0,004x

		double x;
		double y;
		if (currentGeneration < 200) {
			x = 0.004;
			y = 8.9885;
			return (int) Math.round(y * Math.exp(x * -1 * (currentGeneration)));

		} else if (currentGeneration < 700) {
			return 3;
		} else if (currentGeneration < 1000) {
			return 2;
		} else {
			return 1;
		}

	}

	private int getAdaptedLinearTriangleMutation() {
		// f(x) = 8.8929 - 0.031x 
		
		double x;
		double y;
		if (currentGeneration < 200) {
			x = -0.031;
			y = 8.8929;
			return (int) Math.round(currentGeneration * x + y);

		} else if (currentGeneration < 700) {
			return 3;
		} else if (currentGeneration < 1000) {
			return 2;
		} else {
			return 1;
		}
	}

	// Method to select parents using the tournament selection method
	public int[] tournamentSelection() {
		int[] parents = new int[2];
		parents[0] = r.nextInt(populationSize);
		int temp;
		for (int i = 0; i < tournamentSize; i++) {
			temp = r.nextInt(populationSize);
			if (population[temp].getFitness() < population[parents[0]].getFitness()) {
				parents[0] = temp;
			}
		}

		//avoid incest
		parents[1] = r.nextInt(populationSize);
		while (parents[1] == parents[0]) {
			parents[1] = r.nextInt(populationSize);
		}
		for (int i = 0; i < tournamentSize; i++) {
			temp = r.nextInt(populationSize);
			//avoid incest
			while (temp == parents[0]) {
				temp = r.nextInt(populationSize);
			}
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
		
		rouletteWheel.put(0, invertedProbabilities.get(0) / totalFitness);
			
		for (int i = 1; i < populationSize; i++) {
			rouletteWheel.put(i, (invertedProbabilities.get(i) / totalFitness) + rouletteWheel.get(i - 1));
		}

		return rouletteWheel;
	}

	public Solution[] applyCrossover(int[] parents) {
		Solution [] s = new Solution [numberOfOffsprings];
		switch (crossoverMethod) {
		case "standard": return applyStandardCrossover(parents);
		case "PMXO_triangle":
		case "cycle_triangle": 
		case "cycle": 
		case "PMXO": return applyPositionCrossover(parents);
		case "random_triangle": return applyRandomAfterEachTriangleCrossover(parents);
		case "six_way": return applySixWayCrossover(parents);
		case "seperate": return applySeperateCrossover(parents);
		case "optimal": return applyOptimalCrossover (parents);
		case "optimal_mixture": return applyOptimalMixtureCrossover(parents);
		default: return applyStandardCrossover(parents);

		}
	}

	public Solution[] applyStandardCrossover(int[] parents) {
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

	public Solution[] applyPositionCrossover(int[] parents) {
		// initialisation & inversion
		Solution firstParent = population[parents[0]].copy();
		Solution secondParent = population[parents[1]].copy();


		firstParent.applyInversion();
		secondParent.applyInversion();
		Solution offspring[] = new Solution[numberOfOffsprings];

		// crossover
		int windowStartPoint;
		int windowEndPoint;
		int cyclicPointer;
		switch(crossoverMethod) {
			default:
			case "standard": 
				offspring = applyStandardCrossover(parents);
				break;
			case "cycle":
				cyclicPointer = r.nextInt(instance.getNumberOfTriangles()*Solution.VALUES_PER_TRIANGLE);
				offspring [0] = applyCycleCrossover(firstParent, secondParent, cyclicPointer); 
				if (numberOfOffsprings == 2) {
					offspring [1] = applyCycleCrossover(secondParent, firstParent, cyclicPointer);
				}
				break;
			case "PMXO":
				windowStartPoint = r.nextInt(firstParent.getValues().length - 1);
				windowEndPoint = r.nextInt(firstParent.getValues().length-windowStartPoint - 1) + windowStartPoint + 1;
				offspring[0] = applyPartiallyMappedCrossover (firstParent, secondParent, windowStartPoint, windowEndPoint);
				if (numberOfOffsprings == 2) {
					offspring [1] = applyPartiallyMappedCrossover(secondParent, firstParent, windowStartPoint, windowEndPoint);
				}
				break;
			case "cycle_triangle":
				cyclicPointer = r.nextInt(instance.getNumberOfTriangles())*Solution.VALUES_PER_TRIANGLE;
				offspring [0] = applyCycleCrossover(firstParent, secondParent, cyclicPointer); 
				if (numberOfOffsprings == 2) {
					offspring [1] = applyCycleCrossover(secondParent, firstParent, cyclicPointer);
				}
				break;
			case "PMXO_triangle":
				windowStartPoint = r.nextInt(instance.getNumberOfTriangles() - 1);
				windowEndPoint = r.nextInt(instance.getNumberOfTriangles()-windowStartPoint - 1) + windowStartPoint + 1;
				windowStartPoint = windowStartPoint * Solution.VALUES_PER_TRIANGLE;
				windowEndPoint = windowEndPoint * Solution.VALUES_PER_TRIANGLE;
				offspring[0] = applyPartiallyMappedCrossover (firstParent, secondParent, windowStartPoint, windowEndPoint);
				if (numberOfOffsprings == 2) {
					offspring [1] = applyPartiallyMappedCrossover(secondParent, firstParent, windowStartPoint, windowEndPoint);
				}
				break;			
		}

		// reordering
		offspring[0].applyReordering();
		if (numberOfOffsprings == 2) {
			offspring[1].applyReordering();
		}

		return offspring;
	}

	public Solution [] applyRandomAfterEachTriangleCrossover(int[] parents) {
		Solution [] offspring = new Solution [numberOfOffsprings];
		offspring[0] = population[parents[0]].copy();
		if (numberOfOffsprings == 2) offspring [1] = population[parents[1]].copy();
		
		for (int i = 0; i < instance.getNumberOfTriangles() * Solution.VALUES_PER_TRIANGLE; i = i + 10) {
			int currentRandom = r.nextInt(2);
			for (int j = i; j < i + 9; j++) {
				offspring[0].setValue(j, population[parents[currentRandom]].getValue(j));
				if (numberOfOffsprings == 2) {
					offspring[1].setValue(j, population[parents[Math.abs(currentRandom-1)]].getValue(j));
				}
			}
		}
		return offspring;
	}

	public Solution [] applySixWayCrossover(int[] parents) {
		Solution [] offspring = new Solution [numberOfOffsprings];
		offspring [0] = population[parents[0]].copy();
		if (numberOfOffsprings == 2) offspring[1] = population[parents[1]].copy();

		int crossoverPoint1 = r.nextInt(instance.getNumberOfTriangles()) / 4 * Solution.VALUES_PER_TRIANGLE;
		int crossoverPoint2 = (r.nextInt((instance.getNumberOfTriangles() - (crossoverPoint1 / 10)) / 3)
				* Solution.VALUES_PER_TRIANGLE) + crossoverPoint1;
		int crossoverPoint3 = (r.nextInt((instance.getNumberOfTriangles() - (crossoverPoint2 / 10)) / 3)
				* Solution.VALUES_PER_TRIANGLE) + crossoverPoint2;
		int crossoverPoint4 = (r.nextInt((instance.getNumberOfTriangles() - (crossoverPoint3 / 10)) / 3)
				* Solution.VALUES_PER_TRIANGLE) + crossoverPoint3;
		int crossoverPoint5 = (r.nextInt((instance.getNumberOfTriangles() - (crossoverPoint4 / 10)) / 2)
				* Solution.VALUES_PER_TRIANGLE) + crossoverPoint4;
		int crossoverPoint6 = (r.nextInt((instance.getNumberOfTriangles() - (crossoverPoint5 / 10)) / 2)
				* Solution.VALUES_PER_TRIANGLE) + crossoverPoint5;

		int randomParent1 = r.nextInt(2);
		int randomParent2 = r.nextInt(2);
		int randomParent3 = r.nextInt(2);
		int randomParent4 = r.nextInt(2);
		int randomParent5 = r.nextInt(2);
		int randomParent6 = r.nextInt(2);
		int randomParent7 = r.nextInt(2);

		for (int i = 0; i < crossoverPoint1; i++) {
			offspring[0].setValue(i, population[parents[randomParent1]].getValue(i));
			if (numberOfOffsprings == 2) offspring[1].setValue(i, population[parents[Math.abs(randomParent1-1)]].getValue(i));
		}
		for (int i = crossoverPoint1; i < crossoverPoint2; i++) {
			offspring[0].setValue(i, population[parents[randomParent2]].getValue(i));
			if (numberOfOffsprings == 2) offspring[1].setValue(i, population[parents[Math.abs(randomParent2-1)]].getValue(i));

		}
		for (int i = crossoverPoint2; i < crossoverPoint3; i++) {
			offspring[0].setValue(i, population[parents[randomParent3]].getValue(i));
			if (numberOfOffsprings == 2) offspring[1].setValue(i, population[parents[Math.abs(randomParent3-1)]].getValue(i));

		}
		for (int i = crossoverPoint3; i < crossoverPoint4; i++) {
			offspring[0].setValue(i, population[parents[randomParent4]].getValue(i));
			if (numberOfOffsprings == 2) offspring[1].setValue(i, population[parents[Math.abs(randomParent4-1)]].getValue(i));

		}
		for (int i = crossoverPoint4; i < crossoverPoint5; i++) {
			offspring[0].setValue(i, population[parents[randomParent5]].getValue(i));
			if (numberOfOffsprings == 2) offspring[1].setValue(i, population[parents[Math.abs(randomParent5-1)]].getValue(i));

		}
		for (int i = crossoverPoint5; i < crossoverPoint6; i++) {
			offspring[0].setValue(i, population[parents[randomParent6]].getValue(i));
			if (numberOfOffsprings == 2) offspring[1].setValue(i, population[parents[Math.abs(randomParent6-1)]].getValue(i));

		}
		for (int i = crossoverPoint6; i < instance.getNumberOfTriangles() * Solution.VALUES_PER_TRIANGLE; i++) {
			offspring[0].setValue(i, population[parents[randomParent7]].getValue(i));
			if (numberOfOffsprings == 2) offspring[1].setValue(i, population[parents[Math.abs(randomParent7-1)]].getValue(i));

		}
		return offspring;
	}

	public Solution [] applySeperateCrossover(int[] parents) {
		Solution firstParent = population[parents[0]];
		Solution secondParent = population[parents[1]];
		Solution offspring [] = new Solution [numberOfOffsprings];
		offspring [0] = firstParent.copy();
		if (numberOfOffsprings==2) offspring[1] = secondParent.copy();

		for (int selectedtriangle = 0; selectedtriangle < instance.getNumberOfTriangles(); selectedtriangle++) {
			int crossoverpoint = selectedtriangle * 10;
			int crossoverPointColor = crossoverpoint + r.nextInt(4);
			for (int i = crossoverPointColor; i <= crossoverpoint + 3; i++) {
				offspring[0].setValue(i, secondParent.getValue(i));
				if (numberOfOffsprings==2) offspring[1].setValue(i, firstParent.getValue(i));
			}
			int crossoverPointVertex = crossoverpoint + r.nextInt(6) + 4;
			for (int j = crossoverPointVertex; j <= crossoverpoint + 9; j++) {
				offspring[0].setValue(j, secondParent.getValue(j));
				if (numberOfOffsprings==2) offspring[1].setValue(j, firstParent.getValue(j));
			}
		}
		return offspring;
	}

	private Solution[] applyOptimalCrossover(int [] parents) {
		int optimalMethod = 0;
		crossoverMethod = "standard";
		Solution [] optimal = applyCrossover(parents);
		for (int i=0; i<optimal.length; i++) {
			optimal[i].evaluate();
		}
		for (int i=1; i<mutation_methods.length; i++) {
			crossoverMethod = mutation_methods[i];
			Solution [] currentCrossOverSolution = applyCrossover(parents);
			if (improvedFitness(currentCrossOverSolution, optimal)) {
				optimalMethod = i;
				optimal = currentCrossOverSolution;
			}
		}
		 if (currentGeneration < 2000) optimalMethodCount[(int)Math.floor(currentGeneration/optimalCount)][optimalMethod]++;
		crossoverMethod = "optimal";
		return optimal;
	}

	private boolean improvedFitness(Solution[] currentCrossOverSolution, Solution[] optimal) {
		double fitnessCurrent = 0;
		double fitnessOptimal = 0;
		for (int i = 0; i < currentCrossOverSolution.length; i++) {
			currentCrossOverSolution[i].evaluate();
			fitnessCurrent += currentCrossOverSolution[i].getFitness() / currentCrossOverSolution.length;
			fitnessOptimal += optimal[i].getFitness() / currentCrossOverSolution.length;
		}
		if (fitnessCurrent < fitnessOptimal)
			return true;
		return false;
	}
	
	private Solution[] applyOptimalMixtureCrossover(int [] parents) {
		crossoverMethod = "standard";
		Solution [] optimal = applyCrossover(parents);
		for (int i=0; i<optimal.length; i++) {
			optimal[i].evaluate();
		}
		for (int i=1; i<mutation_methods.length; i++) {
			crossoverMethod = mutation_methods[i];
			Solution [] currentCrossOverSolution = applyCrossover(parents);
			optimal = improveFitness(currentCrossOverSolution, optimal);
		}
		crossoverMethod = "optimal_mixture";
		return optimal;
	}

	private Solution [] improveFitness(Solution[] currentCrossOverSolution, Solution[] optimal) {

		for (int i = 0; i<currentCrossOverSolution.length; i++) {
			currentCrossOverSolution[i].evaluate();
			for (int j=0; j<optimal.length; j++) {
				if (currentCrossOverSolution[i].getFitness()<optimal[j].getFitness()) {
					optimal [j] = currentCrossOverSolution[i];
					break;
				}
			}
		}
		return optimal;
	}

	private Solution applyPartiallyMappedCrossover(Solution firstParent, Solution secondParent, int startCrossPoint,
			int endCrossPoint) {
		Solution offspring = firstParent.copy();

		for (int i = 0; i < firstParent.getValues().length; i++) {

			// if value not in window
			if (i < startCrossPoint || i > endCrossPoint) {

				// if orderValue in secondParent window, use first parent's
				// according value
				int indexAtFirstParent = firstParent.getInversionOrderIndex(i);
				if (isInWindow(indexAtFirstParent, secondParent, startCrossPoint, endCrossPoint)) {
					while (isInWindow(indexAtFirstParent, secondParent, startCrossPoint, endCrossPoint)) {
						indexAtFirstParent = firstParent
								.getInversionOrderIndex(secondParent.getInversionArrayPointerIndex(indexAtFirstParent));
					}
					offspring.setValue(i,
							firstParent.getValue(firstParent.getInversionArrayPointerIndex(indexAtFirstParent)));
					// offspring.setInversionOrderIndex(i,
					// firstParent.getInversionOrderIndex(secondParent.getInversionArrayPointerIndex(indexAtFirstParent)));
					offspring.setInversionOrderIndex(i, indexAtFirstParent);

				}

				// else, use first parent's value
				else {
					offspring.setValue(i, firstParent.getValue(i));
					offspring.setInversionOrderIndex(i, indexAtFirstParent);
				}
			}

			// else: copy values from secondParent
			else {
				offspring.setValue(i, secondParent.getValue(i));
				offspring.setInversionOrderIndex(i, secondParent.getInversionOrderIndex(i));
			}
		}

		return offspring;
	}

	private boolean isInWindow(int inversionOrderIndex, Solution secondParent, int startCrossPoint, int endCrossPoint) {
		if (secondParent.getInversionArrayPointerIndex(inversionOrderIndex) < startCrossPoint
				|| secondParent.getInversionArrayPointerIndex(inversionOrderIndex) > endCrossPoint) {
			return false;
		} else {
			return true;
		}
	}

	private Solution applyCycleCrossover(Solution firstParent, Solution secondParent, int cyclicPointer) {

		Solution offspring = firstParent.copy();
		for (int i = 0; i < instance.getNumberOfTriangles() * Solution.VALUES_PER_TRIANGLE; i++) {
			offspring.setValue(i, -1);
		}
		while (offspring.getValue(firstParent.getInversionArrayPointerIndex(cyclicPointer))==-1) {
			cyclicPointer = firstParent.getInversionArrayPointerIndex(cyclicPointer);
			offspring.setValue(cyclicPointer, firstParent.getValue(cyclicPointer));
			offspring.setInversionOrderIndex(cyclicPointer, firstParent.getInversionOrderIndex(cyclicPointer));
			cyclicPointer = secondParent.getInversionOrderIndex(cyclicPointer);
		}

		for (int i = 0; i < instance.getNumberOfTriangles() * Solution.VALUES_PER_TRIANGLE; i++) {
			if (offspring.getValue(i) == -1) {
				offspring.setValue(i, secondParent.getValue(i));
				offspring.setInversionOrderIndex(i, secondParent.getInversionOrderIndex(i));
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
	
	// Calculate diversity measures: Phenotypic, standard deviation
	private void generateDiversityMeasures() {
		double totalFitness = 0.0;
		pVariance = 0.0;
		hammingDistance = 0;
		hammingDistanceRatio = 0.0;
		// Keeps track of the maximum possible distance of the total, which
		maxPossibleDistance = 0;

		for (int i = 0; i < populationSize; i++) {
			totalFitness += population[i].getFitness();
		}
		double averageFitness = totalFitness / populationSize;

		HashSet uniqueFitnesses = new HashSet();
		for (int individual = 0; individual < populationSize; individual++) {
			pVariance += Math.pow(population[individual].getFitness() - averageFitness, 2);
			uniqueFitnesses.add(population[individual].getFitness());

			getHammingDistance(individual);
		}

		hammingDistanceRatio = (double) hammingDistance / maxPossibleDistance;

		pType = uniqueFitnesses.size();

		storeMeasures(averageFitness);

		// Write to file only if last generation of last run has been reached.
		if (Main.currentRun == Main.NUMBER_OF_RUNS - 1 && currentGeneration == Main.NUMBER_OF_GENERATIONS) {
			try {
				writeMeasures();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

	// Stores different measures into arrays so they can be averaged later and
	// saved to a file for analysis.
	private void storeMeasures(double averageFitness) {
		Main.pTypes[Main.currentRun][currentGeneration - 1] = pType;
		Main.hammingDistances[Main.currentRun][currentGeneration - 1] = hammingDistance;
		Main.pVariances[Main.currentRun][currentGeneration - 1] = pVariance;
		Main.hammingDistanceRatios[Main.currentRun][currentGeneration - 1] = hammingDistanceRatio;
		Main.averageFitnesses[Main.currentRun][currentGeneration - 1] = averageFitness;
		Main.bestFitnesses[Main.currentRun][currentGeneration - 1] = getBest(population).getFitness();
	}

	private void getHammingDistance(int individual) {
		for (int secondIndividual = 1; secondIndividual < populationSize; secondIndividual++) {
			if (individual != secondIndividual) {
				for (int valueIndex = 0; valueIndex < Main.NUMBER_OF_TRIANGLES
						* Solution.VALUES_PER_TRIANGLE; valueIndex++) {
					if (population[individual].values[valueIndex] != population[secondIndividual].values[valueIndex]) {
						hammingDistance++;
					}
					maxPossibleDistance++;
				}
			}
		}
	}

	private void writeMeasures() throws IOException {
		// Calculate averages first
		int[] avgPTypeAr = new int[Main.NUMBER_OF_GENERATIONS];
		int[] avgHammingDistanceAr = new int[Main.NUMBER_OF_GENERATIONS];
		double[] avgPVarianceAr = new double[Main.NUMBER_OF_GENERATIONS];
		double[] avgHammingDistanceRatioAr = new double[Main.NUMBER_OF_GENERATIONS];
		double[] avgFitnessAr = new double[Main.NUMBER_OF_GENERATIONS];
		double[] avgBestFitnessAr = new double[Main.NUMBER_OF_GENERATIONS];

		for (int generation = 0; generation < Main.NUMBER_OF_GENERATIONS; generation++) {
			int avgPType = 0;
			int avgHammingDistance = 0;
			double avgPVariance = 0.0;
			double avgHammingDistanceRatio = 0.0;
			double avgFitness = 0.0;
			double avgBestFitness = 0.0;

			for (int run = 0; run < Main.NUMBER_OF_RUNS; run++) {
				avgPType += Main.pTypes[run][generation];
				avgHammingDistance += Main.hammingDistances[run][generation];
				avgPVariance += Main.pVariances[run][generation];
				avgHammingDistanceRatio += Main.hammingDistanceRatios[run][generation];
				avgFitness += Main.averageFitnesses[run][generation];
				avgBestFitness += Main.bestFitnesses[run][generation];
			}
			avgPTypeAr[generation] = avgPType / Main.NUMBER_OF_RUNS;
			avgHammingDistanceAr[generation] = avgHammingDistance / Main.NUMBER_OF_RUNS;
			avgPVarianceAr[generation] = avgPVariance / Main.NUMBER_OF_RUNS;
			avgHammingDistanceRatioAr[generation] = avgHammingDistanceRatio / Main.NUMBER_OF_RUNS;
			avgFitnessAr[generation] = avgFitness / Main.NUMBER_OF_RUNS;
			avgBestFitnessAr[generation] = avgBestFitness / Main.NUMBER_OF_RUNS;

		}
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());

		File file = new File(
				"data/" + timeStamp + " - " + mutationProbability + "-" + numberOfTriangleMutations + ".txt");

		// if file doesnt exists, then create it
		if (!file.exists()) {
			file.createNewFile();
		}

		FileWriter fw = new FileWriter(file.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);

		bw.write(
				"Generation \t PType \t Hamming Distance \t PVariance \t Hamming Distance Ratio \t Avg Fitness \t Best Fitness \n");

		for (int i = 0; i < Main.NUMBER_OF_GENERATIONS; i++) {
			bw.write(i + "\t" + avgPTypeAr[i] + "\t" + avgHammingDistanceAr[i] + "\t"
					+ String.valueOf(avgPVarianceAr[i]).replace('.', ',') + "\t"
					+ String.valueOf(avgHammingDistanceRatioAr[i]).replace('.', ',') + "\t"
					+ String.valueOf(avgFitnessAr[i]).replace('.', ',') + "\t"
					+ String.valueOf(avgBestFitnessAr[i]).replace('.', ',') + "\n");
		}
		bw.close();

	}
}
