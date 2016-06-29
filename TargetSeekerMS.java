import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;

public class TargetSeekerMS {
	private static final String CONTROL_CONDITION = "C";
	private static final String DRUG_CONDITION = "D";
	private static final String REPLICATE_TAG = "R";
	private static final String FRACTION_TAG = "F";
	
	/*
	 * For removing the ProteinID's with the words contaminant and Reverse in
	 * the name e.g. contaminant_GR72_TOBAC Reverse_sp|Q9P281|BAHC1_HUMAN
	 */
	private static final int WINDOW_SIZE = 100; // bins for the distribution
	private static final double STEP_SIZE = .01; // step size for the distribution
	private static final double FPT_STEP_SIZE = .0001;
	private static final int FPT_WINDOW_SIZE = 10000;
	private static final double FPT_SMALLER_STEP_SIZE = .000001;
	// to resolve sticky bit double errors
	private static final double DOUBLE_BUFFER = 0.00000001;
	// the number of conditions should always be 2
	private static final int AVG = 100; 
	// how to signal for Correlation class that you want the average
	private static final int numConditions = 2;
	private static final int SPEC_COUNT_THRESHOLD = 5;
	private static final boolean k_value_DEBUG = false;

	// Stores a list of the ProteinIDs
	private static Set<String> proteinIDSet;
	// a matrix of class Replicate for easy access
	private static Replicate[][] Condition;
	// the new header generated from this program
	private static String[] newHeader;
	// stores and formats the correlation data
	private static Correlation correlationClass;
	
	//stores the distribution, probability, and pvalue items
	private static ResultTable drugTable;
	private static ResultTable[] falsePositiveTable;
	private static double[][] falsePositiveCounter;
	private static double[][] FDRversusPredictions;
	private static Hashtable<String, String> proteinIDtoDescription;
	private static Hashtable<String, Double> foldChange;
	private static ArrayList<String> validProteinIDs;
	private static ArrayList<String> targetProteinIDs;
	
	private static String file_name;				// "Wormdata.txt", "HEKdata.txt", or "Heatdata.txt"
	private static int number_control_replicates;	// 4 or 6 for worm, 4 for HEK, 4 for Heat
	private static int number_drug_replicates;		// 3 for worm, 2 for HEK, 2 for Heat
	private static int number_fractions;			// 10 by default
	private static int k_value;						// 20 by default
	private static double FDR_threshold;			// 0.05 by default
	private static double fold_change_threshold;	// 0.2 by default
	/* controlRepArray and drugRepArray were created to quickly fix an artifact of the
	 * program, involving leaving out certain replicates. (e.g. throwing out replicate 1
	 * out of 6 replicates, would have an array [2,3,4,5,6]) We expect the number of
	 * replicates from our users to be contiguous. Having this properly fixed would not
	 * be worth the time and effort, as it is trivial. For old code dealing with this
	 * artifact, see function setReps().
	 */
	private static String output_file;
	private static int[] controlRepArray;		
	private static int[] drugRepArray;
	private static List<Score> scores;

	

	public TargetSeekerMS(String filePathOnServer, int numControlReplicates, int numDrugReplicates,
					  int numOfFractions, int kValue, double FDRThreshold, double foldChangeThreshold, 
					  String outputFile){
		/* Format: 
		 * 		input_file: .txt file containing spectral counts separated by tabs.
		 * 			File should follow the following naming convention for header:
		 * 			Protein	CR1F1	CR1F2	... DRxFy
		 * 		number_control_replicates: the number of control replicates (integer > 4)
		 * 		number_drug_replicates: the number of drug replicates (integer > 0)
		 * 		number_fractions: the number of fractions in each replicate (integer > 2)
		 * 		k_value: the k-nearest-neighbor smoothing factor (integer > 0)
		 * 		FDR_threshold: the FDR cut-off for determining significant protein.
		 * 			Any proteins with an FDR higher than this threshold are considered
		 * 		 	not significant. 
		 * 		fold_change_threshold: the fold-change threshold for determining 
		 * 			significant protein. Any proteins with a fold-change less than this
		 * 			threshold are considered not significant.
		 */
		
		file_name = filePathOnServer;
		number_control_replicates = numControlReplicates;
		number_drug_replicates = numDrugReplicates;
		number_fractions = numOfFractions;
		k_value = kValue;
		FDR_threshold = FDRThreshold;
		fold_change_threshold = foldChangeThreshold;
		output_file = outputFile;
	}
	
	public void run() {
		
		/* This was how the program was run from the terminal, before it was implemented as an object
		if(args[0].length() < 7){
			System.err.println("Incorrect usage.\nusage: java TargetSeekerMS <input_file>"+
					" <number_control_replicates> <number_drug_replicates> "+
					"<number_fractions> <k_value> <FDR_threshold> <fold_change_threshold>");
			System.exit(-1);
		}
		
		file_name = args[0];
		number_control_replicates = Integer.parseInt(args[1]);
		number_drug_replicates = Integer.parseInt(args[2]);
		number_fractions = Integer.parseInt(args[3]);
		k_value = Integer.parseInt(args[4]);
		FDR_threshold = Double.parseDouble(args[5]);
		fold_change_threshold = Double.parseDouble(args[6]);*/
		
		// set up arrays for control and drug replicates
		controlRepArray = new int[number_control_replicates];
		drugRepArray = new int[number_drug_replicates];
		for(int i = 1; i<=number_control_replicates; i++){
			controlRepArray[i-1] = i;
		}
		for(int i = 1; i<=number_drug_replicates; i++){
			drugRepArray[i-1] = i;
		}


		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(file_name)));
			String s = in.readLine();
			String[] split = null;

			// Stores the file header (starts with Protein)
			String[] header = s.split("\t");

			// Sets up the new file header before spec count
			ArrayList<String> headerColumnsAdded = new ArrayList<String>();
			headerColumnsAdded.add("Protein");
			correlationClass = new Correlation(headerColumnsAdded, controlRepArray, drugRepArray);
			headerColumnsAdded = correlationClass.getHeader();

			// Stores the new, sorted file header.
			newHeader = new String[(headerColumnsAdded.size()
					+ ((controlRepArray.length + drugRepArray.length) * number_fractions))];
			for (int i = 0; i < headerColumnsAdded.size(); i++) {
				newHeader[i] = headerColumnsAdded.get(i);
			}

			// an array of arrays, keeping track of condition and rep number
			Replicate[] Control = new Replicate[controlRepArray.length];
			Replicate[] Drug = new Replicate[drugRepArray.length];
			for (int i = 0; i < controlRepArray.length; i++) { // sets up the
																// control
																// condition
				Control[i] = new Replicate(CONTROL_CONDITION, controlRepArray[i]);
				for (int j = 0; j < number_fractions; j++) {
					String temp = CONTROL_CONDITION + REPLICATE_TAG + controlRepArray[i] +
							FRACTION_TAG + (j + 1);
					newHeader[((number_fractions * i) + j + headerColumnsAdded.size())] = temp;
				}
			}
			for (int i = 0; i < drugRepArray.length; i++) { // sets up the drug
															// condition
				Drug[i] = new Replicate(DRUG_CONDITION, drugRepArray[i]);
				for (int j = 0; j < number_fractions; j++) {
					String temp = DRUG_CONDITION + REPLICATE_TAG + controlRepArray[i] +
							FRACTION_TAG + (j + 1);
					newHeader[((controlRepArray.length * number_fractions) + (number_fractions * i) + j
							+ headerColumnsAdded.size())] = temp;
				}
			}
			Condition = new Replicate[][] { Control, Drug };

			// Finds the order for each rep
			if (header != null) {
				for (int i = 0; i < numConditions; i++) {
					for (int j = 0; j < Condition[i].length; j++) {
						Condition[i][j].setOrder(header, number_fractions);
					}
				}
			}

			// Adds the spec count data
			s = in.readLine();
			while (s != null) {
				split = s.split("\t");

				String ProteinID = split[0];
				
				// i keeps track of the condition number
				// j keeps track of the rep number
				// k keeps track of the frac number
				for (int i = 0; i < numConditions; i++) {
					for (int j = 0; j < Condition[i].length; j++) {
						int[] fracOrder = Condition[i][j].fracOrder();
						List<Double> specData = new ArrayList<Double>();
						for (int k = 0; k < number_fractions; k++) {
							specData.add(Double.parseDouble(split[fracOrder[k]]));
						}

						// loads spec data into the correct rep
						Condition[i][j].loadSpecCounts(ProteinID, specData);
					}
					
				}

				
				s = in.readLine();
			}
			proteinIDSet = Condition[0][0].getProteinIDtoSpecCountVectors().keySet();

			// sets spec and frac abundance
			for (int i = 0; i < controlRepArray.length; i++) { 
				Condition[0][i].setAbundance();
			}
			for (int i = 0; i < drugRepArray.length; i++) { 
				Condition[1][i].setAbundance();
			}
			
			//set spec count threshold
			Condition[0][0].setSpecCountThreshold(SPEC_COUNT_THRESHOLD);
			
			// This normalizes each replicate across the proteinID. This is new.
			normalizeReplicatesAcrossProteinID();
			
			// Computes correlation values
			setUpCorrelationClass();

			computeDrugTable();
			computeFalsePositiveTable();
			drugTable.setFalsePositiveCounter(falsePositiveCounter);
			calculateFoldChange();
			
			computeFDRversusPreictions();
			System.out.println("Done with calculations");

			sortByFDR();
			
			write(output_file);
			in.close();
			System.out.println("done");
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	//unique for this class
	private static void normalizeReplicatesAcrossProteinID() {
		for (int i = 0; i < numConditions; i++) {
			for (int j = 0; j < Condition[i].length; j++) {
				Condition[i][j].normalizeSpecCountAcrossProteinID();
			}
		}
	}
	
	// keeps track of how much the drug changed the protein 
	// flags true if the change is larger than the fold change threshold
	private static void calculateFoldChange() {
		foldChange = new Hashtable<String, Double>();
		
		for(String proteinID: proteinIDSet){
			
			double Corr_cont_avg = correlationClass.getCorrCont(proteinID, AVG, AVG);
			double Corr_D_avg_C_avg = correlationClass.getCorrDCRepXRepY(proteinID, AVG, AVG);
			
			double diff;
			if( (Corr_D_avg_C_avg > 0+DOUBLE_BUFFER) && (Corr_cont_avg > 0-DOUBLE_BUFFER)){
				diff = (Math.abs(Corr_cont_avg-Corr_D_avg_C_avg)/Corr_D_avg_C_avg);
			} else{
				diff = -2;
			}
				
			foldChange.put(proteinID, diff);
				
		}
		drugTable.setFoldChange(foldChange, fold_change_threshold);
		
	}

	// for the purpose of the FDR versus Predictions graph
	// for every FDR value, the number of proteins that have an FDR less than the given value is counted.
	private static void computeFDRversusPreictions() {
		//to choose which FRD method you want. 2 for method 1, 3 for method 2.
		int FDRindex = 2;
		
		// validProteinIDs keeps track of which proteinIDs have valid pValues
		// targetProteinIDs keeps track of the proteinIDs that are under the FDR theshold
		validProteinIDs = new ArrayList<String>();
		targetProteinIDs = new ArrayList<String>();
		
		//a list of unique FDR values, which will be sorted
		ArrayList<Double> FDRList = new ArrayList<Double>();
		FDRList.add(0.0);
		for(double[] d: falsePositiveCounter){
			if(!(FDRList.contains(d[FDRindex]))){
				FDRList.add(d[FDRindex]);
			}
		}
		Collections.sort(FDRList);
		
		int[] FDRCounter = new int[FDRList.size()];
		Hashtable<String, Double> proteinIDtoFDR = new Hashtable<String, Double>();

		for(String proteinID: proteinIDSet){
			//pValue is necessary for grabbing the FDR
			double pValue = drugTable.getProteinNoisePValue().get(proteinID);
			double FDRValue = 0;

			if(pValue == -2.0){
				proteinIDtoFDR.put(proteinID, -2.0);
			}
			if( (pValue > (0-DOUBLE_BUFFER)) && (pValue < DOUBLE_BUFFER) ){
				//pValues are listed in decreasing order
				proteinIDtoFDR.put(proteinID,falsePositiveCounter[falsePositiveCounter.length-1][FDRindex]);
			} else{
				for(int i = 0; i < falsePositiveCounter.length-1; i++){
					// this rounds them up. This is okay because you should not have a 0 pValue
					if( (pValue < falsePositiveCounter[i][0]+DOUBLE_BUFFER) && (pValue > falsePositiveCounter[i+1][0]+DOUBLE_BUFFER) ){
						//prints everything in the FPC but the column with the PValue bin
						proteinIDtoFDR.put(proteinID,falsePositiveCounter[i][FDRindex]);
						FDRValue = falsePositiveCounter[i][FDRindex];
					} 
				}	
			}
			
			
			if( (pValue >= -(0+DOUBLE_BUFFER)) && (pValue <= (1.0+DOUBLE_BUFFER))){ 
				// adds the protein to the list of valid proteinIDs if it has a valid pValue
				validProteinIDs.add(proteinID);
				
				//WE do not know what this does
				if( (pValue > (0-DOUBLE_BUFFER)) && (pValue < DOUBLE_BUFFER) ){
					// special case for FDR of 0
					FDRValue = 0;
				}
				
				// iterates every value in the FDR counter lower than the current FDR value
				for(int i = FDRList.size()-1; i >= 0; i--){
					if( FDRValue <= FDRList.get(i) ){
						//System.out.println("FDRValue = "+FDRValue+"\nFDRList.get(i) = "+FDRList.get(i));
						FDRCounter[i]++;
					} else{
						break;
					}
				}
				
			} 
			
			if(proteinIDtoFDR.get(proteinID)==null){
				System.out.println(pValue);
			}
			double proteinFDR = proteinIDtoFDR.get(proteinID);
			//adds to the target protein list
			if( (proteinFDR < FDR_threshold) && (proteinFDR > (0-DOUBLE_BUFFER)) &&
						(foldChange.get(proteinID)>=fold_change_threshold)){
					targetProteinIDs.add(proteinID);
			}
		}
		
		// puts the working values into the correct format
		FDRversusPredictions = new double[FDRList.size()][];
		for(int i = 0; i < FDRList.size(); i++){
			FDRversusPredictions[i] = new double[]{FDRList.get(i), FDRCounter[i]};
		}
	}

	// determines how many control and how many drug reps there are
/*	private static void setReps(String[] header) {
		if(PRESET_REPS){ // pre-determined number of reps
			controlRepArray = new int[]{2,3,4,5};
			drugRepArray = new int[]{2,3}; 
			
			// makes the description accessible
			descriptionIndex = 0;
			for (int i = 0; i < header.length; i++) {
				String s = header[i];
				if(s.contains(DESCRIPTION)){
					descriptionIndex = i;
				}
			}

		} else{ // determines how many reps from the file header
			List<Integer> controlReplicates = new ArrayList<Integer>();
			List<Integer> drugReplicates = new ArrayList<Integer>();
			descriptionIndex = 0;

			for (int i = 0; i < header.length; i++) {
				String s = header[i];
				// only considered if it contains "spectral_count"
				if (s.contains(SPEC_TAG)) {
					// divides the spec header by _
					// will look like: spectral_count_Cond_X_Rep_Y_Frac_Z"
					// 0 1 2 3 4 5 6 7
					String[] dividedSpecHeader = s.split("_");
					String condition = dividedSpecHeader[3];
					String replicate = dividedSpecHeader[5];

					if (condition.equals(CONTROL_CONDITION)) {
						if (!(controlReplicates.contains(Integer.parseInt(replicate)))) {
							controlReplicates.add(Integer.parseInt(replicate));
						}
					} else if (condition.equals(DRUG_CONDITION)) {
						if (!(drugReplicates.contains(Integer.parseInt(replicate)))) {
							drugReplicates.add(Integer.parseInt(replicate));
						}
					}
				}
				
				if(s.contains(DESCRIPTION)){
					descriptionIndex = i;
				}
			}

			// makes sure they are ascending order
			Collections.sort(controlReplicates);
			Collections.sort(drugReplicates);

			// adds elements into the array
			for (int i = 0; i < controlReplicates.size(); i++) {
				controlRepArray[i] = controlReplicates.get(i);
			}

			for (int i = 0; i < drugReplicates.size(); i++) {
				drugRepArray[i] = drugReplicates.get(i);
			}
		}
	}*/

	// All correlation calculations are done here
	private static void setUpCorrelationClass() {
		correlationClass.initializeSet(proteinIDSet);

		for (String proteinID : proteinIDSet) {

			/*
			 * Computes the correlation values between drug and control.
			 * Annotated as Corr_DC_Rep_x where x is the rep number
			 */
			if (correlationClass.isCorrDCRep()) {
				for (int j = 0; j < Math.min(controlRepArray.length, drugRepArray.length); j++) {
					List<Double> specCountsCond1 = Condition[0][j].getProteinIDtoSpecCountVectorsNormalizedAcrossProteinID().get(proteinID);
					List<Double> specCountsCond2 = Condition[1][j].getProteinIDtoSpecCountVectorsNormalizedAcrossProteinID().get(proteinID);
					double corrVal = computeCorrelation(specCountsCond1, specCountsCond2);
					correlationClass.addCorrDCRep(proteinID, corrVal, controlRepArray[j]);
				}
			}

			/*
			 * Computes the correlation values between the control reps and its
			 * average Annotated as Corr_cont_x_y where x and y are the two reps
			 * acted on Annotated as Corr_cont_avg for the average If
			 * applicable, the average excluding a Corr_cont_x_y is computed
			 * Annotated as Corr_cont_avg_EXCLUDE_x_y
			 */
			if (correlationClass.isCorrContXY()) {
				Hashtable<String, Double> correlationTable = new Hashtable<String, Double>();
				for (int i = 0; i < controlRepArray.length; i++) {
					for (int j = (i + 1); j < controlRepArray.length; j++) {
						List<Double> specCountsContRepX = Condition[0][i].getProteinIDtoSpecCountVectorsNormalizedAcrossProteinID().get(proteinID);
						List<Double> specCountsContRepY = Condition[0][j].getProteinIDtoSpecCountVectorsNormalizedAcrossProteinID().get(proteinID);
						double corrContX_Y = computeCorrelation(specCountsContRepX, specCountsContRepY);
						String tempName = ("Corr_cont_" + controlRepArray[i] + "_" + controlRepArray[j]);
						correlationTable.put(tempName, corrContX_Y);
					}
				}
				correlationClass.addCorrCont(proteinID, correlationTable);
			}

			/*
			 * Computes the correlation between the condition 2 rep 1 and the
			 * control reps Annotated as Corr_D_Rep1_C_Repx where x is the
			 * control rep number Annotated as Corr_D_Rep1_C_avg for the average
			 */
			if (correlationClass.isCorrDRepCRep()) {
				Hashtable<String, Double> correlationTable = new Hashtable<String, Double>();
				for (int i = 0; i < drugRepArray.length; i++) {
					for (int j = 0; j < controlRepArray.length; j++) {
						List<Double> specCountsDrugRepX = Condition[1][i].getProteinIDtoSpecCountVectorsNormalizedAcrossProteinID().get(proteinID);
						List<Double> specCountsContRepY = Condition[0][j].getProteinIDtoSpecCountVectorsNormalizedAcrossProteinID().get(proteinID);
						String tempName = ("Corr_DrugRep_" + drugRepArray[i] + "_ContRep_" + controlRepArray[j]);
						double corrDX_CY = computeCorrelation(specCountsDrugRepX, specCountsContRepY); // Corr_DrugRepX_ContRepY
						correlationTable.put(tempName, corrDX_CY);
					}
				}
				correlationClass.addCorrDRepCRep(proteinID, correlationTable);
			}

		}

	}

	// adds frequency to the correct bin
	private static double[] addFrequency(double[] distribution, double correlation) {
		// Three way split:
		// ** From -1.00 to (but not including) -STEP_SIZE
		// ** From -STEP_SIZE to STEP_SIZE
		// ** From STEP_SIZE (but not including) to 1.00
		// Steps through each bin and increments if the correlation is in the
		// bin. For WINDOW_SIZE 100 and STEP_SIZE .1 :
		// **bins 0 to 99, bin 100, bins 101 to 200. 201 bins total
		
		//CHANGED 2-16-16 only from 0 to 1.00 now
		// ** bins 0 to 99. 100 total

		if ((correlation>=(0.0-DOUBLE_BUFFER))&&(correlation<=(1.0+DOUBLE_BUFFER))){
			for (int i = 0; i < WINDOW_SIZE; i++) {
				double lowerBound;
				double upperBound;
				
				if(i == 0){
					lowerBound = (0.0-DOUBLE_BUFFER);
				} else {
					lowerBound = 0.00 + (STEP_SIZE * (i));
				}
				
				if( i == WINDOW_SIZE-1){
					upperBound = (1.0+DOUBLE_BUFFER);
				} else {
					upperBound = 0.00 + (STEP_SIZE * (i+1) );	
				}
				
				if ((correlation > lowerBound) && (correlation <= upperBound)) {
					distribution[i]++;
					//break; 
					//check if break is okay. I guess it isn't.
				}
			}
		} 
		return distribution;
	}

	private static void computeDrugTable() {
		/*
		 * Computing frequency matrix, smoothed frequency matrix,
		 * probability matrix average distribution, average probability,
		 * conditional distribution, conditional probability, pvalue
		 */
		drugTable = new ResultTable(true);
		
		//set up correlation tables
		Hashtable<String, double[]> excludeTable = new Hashtable<String, double[]>();
		Hashtable<String, double[]> includeTable = new Hashtable<String, double[]>();
		Hashtable<String, Double> targetTable = new Hashtable<String, Double>();
		int size = ((controlRepArray.length)*(controlRepArray.length-1)/2);
		for(String proteinID: proteinIDSet){
			double[] exclude = new double[size];
			double[] include = new double[size];
			// keeps track of where to add the value
			int count = 0; 
			for(int i = 0; i < controlRepArray.length; i++){
				for(int j = (i+1); j < controlRepArray.length; j++){
					exclude[count] = correlationClass.getCorrContExclude(proteinID, controlRepArray[i], controlRepArray[j]);
					include[count] = correlationClass.getCorrCont(proteinID, controlRepArray[i], controlRepArray[j]);
					count++;
				}
			}
			excludeTable.put(proteinID, exclude);
			includeTable.put(proteinID, include);
			
			//adds target value
			targetTable.put(proteinID, correlationClass.getCorrDCRepXRepY(proteinID, AVG, AVG));
		}
		
		drugTable.setExcludeTable(excludeTable);
		drugTable.setIncludeTable(includeTable);
		drugTable.setTargetTable(targetTable);
		drugTable.setFrequencyMatrix(computeFrequencyMatrix(excludeTable, includeTable));
		drugTable.setSmoothedFrequencyMatrix(smoothDistribution(drugTable.getFrequencyMatrix(), k_value));
		drugTable.setProbabilityMatrix(computeProbabilityMatrix(drugTable.getSmoothedFrequencyMatrix()));
		drugTable.setAverageDistribution(computeAverageDistribution(drugTable.getSmoothedFrequencyMatrix()));
		drugTable.setAverageProbability(computeProbability(drugTable.getAverageDistribution()));
		drugTable.setProteinNoiseProbability(computeProteinNoiseProbability(
				drugTable.getAverageProbability(), drugTable.getProbabilityMatrix()));
		drugTable.setProteinNoisePValue(computePValue(drugTable.getProteinNoiseProbability(), targetTable));
		drugTable.setProteinIDtoDescription(proteinIDtoDescription);
		
		
	}

	private static void computeFalsePositiveTable() {
		
		falsePositiveTable = new ResultTable[controlRepArray.length];
		
		for(int n = 0; n < falsePositiveTable.length; n++){
			// controlRepArray[n] is the control rep we are using as the "drug"
			ResultTable falsePositive = new ResultTable(false);
			
			//set up correlation tables
			Hashtable<String, double[]> excludeTable = new Hashtable<String, double[]>();
			Hashtable<String, double[]> includeTable = new Hashtable<String, double[]>();
			Hashtable<String, Double> targetTable = new Hashtable<String, Double>();
			int size = ((controlRepArray.length-1)*(controlRepArray.length-2)/2);
			for(String proteinID: proteinIDSet){
				double[] exclude = new double[size];
				double[] include = new double[size];
				// keeps track of where to add the value
				int count = 0; 
				for(int i = 0; i < controlRepArray.length; i++){
					for(int j = (i+1); j < controlRepArray.length; j++){
						//only if either i or j is not the same as the control rep we're excluding
						if((controlRepArray[i] != controlRepArray[n]) && (controlRepArray[j] != controlRepArray[n]) ){
							exclude[count] = correlationClass.getCorrContExclude(proteinID, controlRepArray[i], controlRepArray[j]);
							include[count] = correlationClass.getCorrCont(proteinID, controlRepArray[i], controlRepArray[j]);
							count++;
						}
					}
				}
				excludeTable.put(proteinID, exclude);
				includeTable.put(proteinID, include);
				
				//adds target value
				double sum = 0.0;
				int numElements = 0;
				for(int i = 0; i < controlRepArray.length; i++){
					//access corr_cont_x_y values for everything related to "n"
					if( i > n ){
						double value = correlationClass.getCorrCont(proteinID, controlRepArray[n], controlRepArray[i]);
						if( (value > -(1+DOUBLE_BUFFER)) && (value < (1+DOUBLE_BUFFER))){
							sum += value;
							numElements++;
						}
					}
					if( i < n ){
						double value = correlationClass.getCorrCont(proteinID, controlRepArray[i], controlRepArray[n]);
						if( (value > -(1+DOUBLE_BUFFER)) && (value < (1+DOUBLE_BUFFER))){
							sum += value;
							numElements++;
						}
					} 
				}
				targetTable.put(proteinID, (sum/numElements));
				
			}
	
			falsePositive.setExcludeTable(excludeTable);
			falsePositive.setIncludeTable(includeTable);
			falsePositive.setTargetTable(targetTable);
			falsePositive.setFrequencyMatrix(computeFrequencyMatrix(excludeTable, includeTable));
			falsePositive.setSmoothedFrequencyMatrix(smoothDistribution(falsePositive.getFrequencyMatrix(), k_value));
			falsePositive.setProbabilityMatrix(computeProbabilityMatrix(falsePositive.getSmoothedFrequencyMatrix()));
			falsePositive.setAverageDistribution(computeAverageDistribution(falsePositive.getSmoothedFrequencyMatrix()));
			falsePositive.setAverageProbability(computeProbability(falsePositive.getAverageDistribution()));
			falsePositive.setProteinNoiseProbability(computeProteinNoiseProbability(
					falsePositive.getAverageProbability(), falsePositive.getProbabilityMatrix()));
			falsePositive.setProteinNoisePValue(computePValue(falsePositive.getProteinNoiseProbability(), targetTable));
			falsePositiveTable[n] = falsePositive;
		}
		
		//sets up the false positive counter 
		int fptSize = ((FPT_WINDOW_SIZE)+(int)(FPT_STEP_SIZE/FPT_SMALLER_STEP_SIZE));
		falsePositiveCounter = new double[fptSize][];
		for(int i = 0; i < falsePositiveCounter.length; i++){
			if(i < FPT_WINDOW_SIZE){
				falsePositiveCounter[i] = new double[4];// magic number
				falsePositiveCounter[i][0] = (1-(i*FPT_STEP_SIZE));
			}else{
				falsePositiveCounter[i] = new double[4];// magic number
				falsePositiveCounter[i][0] = (FPT_STEP_SIZE - ((i+1)-FPT_WINDOW_SIZE)*FPT_SMALLER_STEP_SIZE);
			}
		}

		//Counts the false positives
		for(int n = 0; n< falsePositiveTable.length; n++){
			Hashtable<String, Double> pValueSet = falsePositiveTable[n].getProteinNoisePValue();
			
			//for method 2
			double[][] individualFalsePositiveCounter = new double[falsePositiveCounter.length][];
			//copies the pValue size of the original array and zeroes out the counter part
			for(int i = 0; i < falsePositiveCounter.length; i++){
				individualFalsePositiveCounter[i] = Arrays.copyOf(falsePositiveCounter[i], 2);
				individualFalsePositiveCounter[i][1] = 0;
			}
			
			for(String proteinID: proteinIDSet){
				Double pValue = pValueSet.get(proteinID);
				
				//excludes invalid pValues
				if((pValue >= (0-DOUBLE_BUFFER)) && (pValue <= (1+DOUBLE_BUFFER) )){
					//adds the pValue to all applicable bins
					for(int i = 0; i < falsePositiveCounter.length-1; i++){
						if(pValue < falsePositiveCounter[i][0]+DOUBLE_BUFFER){
							falsePositiveCounter[i][1]++;
							individualFalsePositiveCounter[i][1]++;
						}else{
							break;
						}
					}
				}
			}
			
			//for method 2
			falsePositiveTable[n].setFalsePositiveCounter(individualFalsePositiveCounter);
		}
		
		//for methods 1 and 2, the denominator
		double[][] truePositiveCounter = new double[falsePositiveCounter.length][];
		//copies the pValue size of the original array and zeroes out the counter part
		for(int i = 0; i < truePositiveCounter.length; i++){
			truePositiveCounter[i] = Arrays.copyOf(falsePositiveCounter[i], 2);
			truePositiveCounter[i][1] = 0;
		}
		for(String proteinID: proteinIDSet){
			Double pValue = drugTable.getProteinNoisePValue().get(proteinID);
			//excludes invalid pValues
			if((pValue >= (0-DOUBLE_BUFFER) ) && (pValue <= (1+DOUBLE_BUFFER))){
				//adds the pValue to all applicable bins
				for(int i = 0; i < truePositiveCounter.length; i++){
					if(pValue <= truePositiveCounter[i][0]+DOUBLE_BUFFER){
						truePositiveCounter[i][1]++;
					}else{
						break;
					}
				}
			}
		}

		//for debugging
		int[] debugQueue = new int[]{0, 1, 50, 500, 5000};
		if(k_value_DEBUG){
			System.out.println("k = "+k_value);
		}
		
		//method 1
		int totalValidPValues = 0;
		double method1MinValue = Double.MAX_VALUE; // for making the function monotonic
		for(int n = 0; n < falsePositiveTable.length; n++){
			totalValidPValues += falsePositiveTable[n].getNumValidPValues();
		}
		
		if(k_value_DEBUG){
			System.out.println("totalValidPValues = "+totalValidPValues);
		}
		
		for(int i = 0; i < falsePositiveCounter.length; i++){
			double numerator = (falsePositiveCounter[i][1]/totalValidPValues);
			double denominator = (truePositiveCounter[i][1]/drugTable.getNumValidPValues());
			double result;

			if(k_value_DEBUG){				
				if( i == 0 ){
					for(ResultTable a: falsePositiveTable){
						System.out.println("FPC[0][1] = "+a.getFalsePositiveCounter()[0][1]+"\tnumValidPValues = "+a.getNumValidPValues());
					}
				}
				for(int d: debugQueue){
					if( i == d )
						System.out.println("For i = "+i+" method 1, numerator: "+numerator+" ; denominator : "+denominator);
				}
			}
			
			//to prevent result is NaN or infinity
			if( (numerator < (0+DOUBLE_BUFFER)) || (denominator < (0)+DOUBLE_BUFFER)){
				result = 0;
			}else{
				result = numerator/denominator;
			}
			
			//monotonic transformation
			if( result < method1MinValue ){
				method1MinValue = result;
			}
			falsePositiveCounter[i][2] = method1MinValue;
		}
		
		//method 2 
		//averages the ratios of false positives to total valid p values
		int[] validPValues = new int[falsePositiveTable.length];
		double method2MinValue = Double.MAX_VALUE; // for making the function monotonic
		for(int n = 0; n < falsePositiveTable.length; n++){
			validPValues[n] = falsePositiveTable[n].getNumValidPValues();
		}


		if(k_value_DEBUG){
			System.out.println("validPValues: ");
			for(int d: validPValues){
				System.out.print(d+"\t");
			}
			System.out.println();
		}
		
		for(int i = 0; i< falsePositiveCounter.length; i++){
			double sum = 0;
			int numElements = 0;
			for(int n = 0; n < falsePositiveTable.length; n++){
				double ratio = (falsePositiveTable[n].getFalsePositiveCounter()[i][1])/(validPValues[n]);
				sum += ratio;
				numElements++;
			}
			
			double numerator = (sum/numElements);
			double denominator = (truePositiveCounter[i][1]/drugTable.getNumValidPValues());
			double result;

			if(k_value_DEBUG){
				for(int d: debugQueue){
					if( i == d )
						System.out.println("For i = "+i+" method 2, numerator: "+numerator+" ; denominator : "+denominator);

				}
			}
			
			//to prevent result is NaN or infinity
			if( (numerator < (0+DOUBLE_BUFFER)) || (denominator < (0)+DOUBLE_BUFFER)){
				result = 0;
			}else{
				result = numerator/denominator;
			}
			
			// for monotonicity
			if( result < method2MinValue ){
				method2MinValue = result;
			}
			falsePositiveCounter[i][3] = method2MinValue;
		}
		
		
	}

	// computes the frequency matrix
	private static double[][] computeFrequencyMatrix(Hashtable<String, double[]> excludeTable,
				Hashtable<String, double[]> includeTable) {
		double[][] frequencyMatrix = new double[WINDOW_SIZE][WINDOW_SIZE];
		
		for (String proteinID : proteinIDSet) { 
			//exclude table holds all the values of Corr_cont_avg_EXCLUDE_x_y
			//include table holds all the values of Corr_cont_x_y
			//the x and y values for these must line up
			double[] exclude = excludeTable.get(proteinID);
			double[] include = includeTable.get(proteinID);
			//exclude stands for the average correlation (excluding the value of include)
			//include is the correlation value we are adding for a given average (exclude)

			for (int i = 0; i < exclude.length; i++) {
				// finds the bin for exclude
				if ((exclude[i]>=(0.0-DOUBLE_BUFFER))&&(exclude[i]<=(1.0+DOUBLE_BUFFER))){
					for (int k = 0; k < WINDOW_SIZE; k++) {
						double lowerBound;
						double upperBound;
						
						if(k == 0){
							lowerBound = (0.0-DOUBLE_BUFFER);
						} else {
							lowerBound = 0.00 + (STEP_SIZE * (k));
						}
						
						if( k == WINDOW_SIZE-1){
							upperBound = (1.0+DOUBLE_BUFFER);
						} else {
							upperBound = 0.00 + (STEP_SIZE * (k+1));	
						}
						
						//System.out.println("lowerBound = "+lowerBound+" \tupperBound = "+upperBound);
						
						if ((exclude[i] > lowerBound) && (exclude[i] <= upperBound)) {
							frequencyMatrix[k] = addFrequency(frequencyMatrix[k], include[i]); // adds
																							// include
						}			
						
					}
				} 
			}
		}
		//it's y and x because thought it appears as x and y on the frequency matrix,
		//the indices are actually swapped on the data structure
		//System.out.println("frequencyMatrix["+y+"]["+x+"] = "+frequencyMatrix[y][x]);
		
		return frequencyMatrix;
	}

	// Smoothing method from Mathieu
	private static double[][] smoothDistribution(double[][] distribution, int k) {
		// Creates new array that will contain the smoothed frequencies
		double[][] smoothedDistribution = new double[distribution.length][distribution.length];
		double[][] sums = new double[distribution.length][distribution.length];
		for (int i = 0; i < sums.length; i++) {
			for (int j = 0; j < sums.length; j++) {
				if (i == 0 && j == 0) {
					sums[i][j] = distribution[i][j];
				} else if (i == 0) {
					sums[i][j] = distribution[i][j] + sums[i][j - 1];
				} else if (j == 0) {
					sums[i][j] = distribution[i][j] + sums[i - 1][j];
				} else {
					sums[i][j] = distribution[i][j] + sums[i][j - 1] + sums[i - 1][j] - sums[i - 1][j - 1];
				}
			}
		}
		for (int i = 0; i < distribution.length; i++) {
			for (int j = 0; j < distribution.length; j++) {
				boolean hasK = false; // Flag that shows if k datapoints have
										// been seen
				if (distribution[i][j] >= k) {
					smoothedDistribution[i][j] = distribution[i][j];
					hasK = true;
				}
				int iUL = i;
				int iUR = i;
				int iLL = i;
				int iLR = i;
				int jUL = j;
				int jUR = j;
				int jLL = j;
				int jLR = j;
				double numberOfCells = 1;

				int prevCount = (int) distribution[i][j]; // gives the amount of
															// frequency seen
															// already
				int prevNumberOfCells = 1;
				int count = 0;
				while (hasK == false) {
					if (iUL > 0) {
						iUL--;
					}
					if (iUR > 0) {
						iUR--;
					}
					if (iLL < distribution.length - 1) {
						iLL++;
					}
					if (iLR < distribution.length - 1) {
						iLR++;
					}
					if (jUL > 0) {
						jUL--;
					}
					if (jUR < distribution.length - 1) {
						jUR++;
					}
					if (jLL > 0) {
						jLL--;
					}
					if (jLR < distribution.length - 1) {
						jLR++;
					}

					if (jLL - 1 < 0) {
						// TOP LEFT
						if (iUR - 1 < 0) {
							count = (int) sums[iLR][jLR];
						}
						// LEFT
						else {
							count = (int) (sums[iLR][jLR] - sums[iUR - 1][jUR]);
						}
					} else if (iUR - 1 < 0) {
						count = (int) (sums[iLR][jLR] - sums[iLL][jLL - 1]);
					} else {
						count = (int) (sums[iLR][jLR] - sums[iLL][jLL - 1] - sums[iUR - 1][jUR]
								+ sums[iUL - 1][jUL - 1]);
					}
					numberOfCells = (jUR - jUL + 1) * (iLL - iUL + 1);

					if (count == k) {
						double weights = (double) 1.0 / ((double) (1.0 * numberOfCells));
						// weighs values by how far they are from the original
						// point in the 2d array
						smoothedDistribution[i][j] = weights * count;
						hasK = true;
					} else if (count > k) {
						double borderCount = count - prevCount;
						double borderNumberCells = numberOfCells - prevNumberOfCells;

						double w2 = (double) (k - prevCount) / (double) (borderCount);
						double weights = 1.0 / ((double) (w2 * borderNumberCells) + (double) (1 * prevNumberOfCells));
						double sum = prevCount + w2 * borderCount;
						smoothedDistribution[i][j] = weights * sum;
						hasK = true;
					}
					prevCount = count;
					prevNumberOfCells = (int) numberOfCells;
				}
			}
		}

		return smoothedDistribution;
	}

	// probability matrix is computed by computing the probability for each row
	// of the smoothed frequency matrix
	private static double[][] computeProbabilityMatrix(double[][] smoothedFrequencyMatrix) {
		double[][] probabilityMatrix = new double[smoothedFrequencyMatrix.length][smoothedFrequencyMatrix[0].length];
		for (int i = 0; i < probabilityMatrix.length; i++) {
			probabilityMatrix[i] = computeProbability(smoothedFrequencyMatrix[i]);
		}
		return probabilityMatrix;
	}

	// collapses the frequency matrix for each row
	private static double[] computeAverageDistribution(double[][] frequencyMatrix) {
		double[] distribution = new double[WINDOW_SIZE];
		for (int i = 0; i < frequencyMatrix.length; i++) {
			double sum = 0;
			for (int j = 0; j < frequencyMatrix.length; j++) {
				sum += frequencyMatrix[i][j];
			}
			distribution[i] = sum;
		}
		return distribution;
	}

	// normalizes the distribution
	private static double[] computeProbability(double[] distribution) {
		double[] probability = new double[distribution.length];
		double sum = 0;
		// finds sum of distribution
		for (double count : distribution) {
			sum += count;
		}
		// calculates probability
		if (sum != 0) {
			for (int i = 0; i < distribution.length; i++) {
				probability[i] = distribution[i] / sum;
			}
		}
		return probability;
	}

	// determines the probability that given a target probability, that the
	// change is due to noise
	private static Hashtable<String, double[]> computeProteinNoiseProbability(double[] averageProbability,
			double[][] probabilityMatrix) {
		Hashtable<String, double[]> proteinNoiseProbability = new Hashtable<String, double[]>();

		for (String proteinID : proteinIDSet) {
			double[] conditionalProbability = new double[WINDOW_SIZE];
			double[] conditionalProbabilityNormalized = new double[WINDOW_SIZE];
			// stores bin number for each corr_cont_x_y
			ArrayList<Integer> binList = new ArrayList<Integer>();

			for (int i = 0; i < controlRepArray.length; i++) {
				for (int j = (i + 1); j < controlRepArray.length; j++) {
					binList.add(
						getBinNum(correlationClass.getCorrCont(proteinID, controlRepArray[i], controlRepArray[j])));
				}
			}

			double sum = 0;
			boolean isNA = false;

			for (int i : binList) { // if any value is -2, it is not valid
				if (i < 0 )
					isNA = true;
			}

			// steps through each bin to compute distribution
			for (int i = 0; i < conditionalProbability.length; i++) {
				if (!isNA) {
					double probability = averageProbability[i];
					for (int j : binList) {
						probability *= probabilityMatrix[i][j];
					}
				    
					conditionalProbability[i] = probability;
				} else {
					conditionalProbability[i] = -2.00;
				}
				sum += conditionalProbability[i];
			}
			if (!isNA) { // normalizes each distribution
				for (int i = 0; i < conditionalProbability.length; i++) {
					conditionalProbabilityNormalized[i] = conditionalProbability[i] / sum;
				}
			} else {
				for (int i = 0; i < conditionalProbability.length; i++) {
					conditionalProbabilityNormalized[i] = -2.00; // to print N/A
																	// later
				}
			}

			proteinNoiseProbability.put(proteinID, conditionalProbabilityNormalized);
		}

		return proteinNoiseProbability;
	}

	// sums up all the values up to a target value for the noise value
	private static Hashtable<String, Double> computePValue(Hashtable<String, double[]> proteinNoiseValue,
				Hashtable<String, Double> targetTable) {
		Hashtable<String, Double> output = new Hashtable<String, Double>();

		for (String proteinID : proteinIDSet) {
			double[] range = proteinNoiseValue.get(proteinID); // distribution
																// of noise
																// values
			double targetValue = targetTable.get(proteinID);
			double pValue = 0;

			if ((range[0] < (0.00-DOUBLE_BUFFER)) || (targetValue < -(0.00+DOUBLE_BUFFER))) {
				pValue = -2.0;
			}else{
				for (int i = 0; i <= getBinNum(targetValue); i++) {
					pValue += range[i];
				}	
			}
			

			output.put(proteinID, pValue);
			
		}

		return output;
	}

	private static int getBinNum(double correlation) {

		if (correlation != -2.0) { // if not -2.0 then write computed
									// correlation
			if(correlation == 1.0){
			  return 99;
			}
			
			int bin = (int) (Math.ceil((correlation) * 100)-1);
			if(bin == -1){
				return 0;
			}
			else{
				return bin;
			}
			
		} else { // if -2.0 then write N/A
			return -2;
		}

	}

	private static double computeCorrelation(List<Double> specCountsCond1, List<Double> specCountsCond2) { 
		//unique to this class. Actually computes 1-Euclidean distance
		double euclidianDistance = 0;
		double sum = 0;
		boolean hasNA = false;
		
		for(int i = 0; i < specCountsCond1.size(); i++){
			if((specCountsCond1.get(i) > -(0+DOUBLE_BUFFER)) && (specCountsCond2.get(i) > -(0+DOUBLE_BUFFER))){
				double diff = (specCountsCond1.get(i)-specCountsCond2.get(i));
				sum += (diff*diff);
			} else{
				hasNA = true;
			}
		}
		euclidianDistance = Math.sqrt(sum);
		
		//normalizes by max distance
		//makes the distance metric between 0 and 1
		euclidianDistance = euclidianDistance/(Math.sqrt(2));
		
		if(hasNA){
			return -2.0;
		} else{
			//transforms the distance to a similarity measure
			return (1.0-euclidianDistance);
		}
	}

	// Sorts by FDR
	private static void sortByFDR(){
		scores = new ArrayList<Score>();
		Hashtable<String,Double> pValueTable = drugTable.getProteinNoisePValue();
			
		for(String proteinID: proteinIDSet){	
			//PValue and FDR
			double pValue = pValueTable.get(proteinID);
			//FDR is 2, so it is sorted accordingly
			double FDRvalue = 2;
			
			if( (pValue > (0-DOUBLE_BUFFER)) && (pValue <= (1.0+DOUBLE_BUFFER))){
				if( (pValue > (0-DOUBLE_BUFFER)) && (pValue < DOUBLE_BUFFER) ){
					FDRvalue = falsePositiveCounter[falsePositiveCounter.length-1][2];
				} else{
					for(int i = 0; i < falsePositiveCounter.length-1; i++){
						if( (pValue < falsePositiveCounter[i][0]+DOUBLE_BUFFER) && (pValue > falsePositiveCounter[i+1][0]+DOUBLE_BUFFER) ){
							FDRvalue = falsePositiveCounter[i][2];
						}
					}	
				}
						
			}
			scores.add(new Score(FDRvalue,proteinID));
		}
		Collections.sort(scores);
	
	}

	/*
	 * Writes the output file. You can toggle whether to print: - spec values -
	 * frequency matrix - smoothed frequency matrix - probability matrix -
	 * average distribution - average probability - rejected proteins
	 */
	private static void write(String outputFile) {
		System.out.println("in write");
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
	
			out.write("Protein\tSimilarity_Drug_vs_Control\tpValue\tFDR\t");
			out.write("Fold-change\tSignificant\n");
			out.flush();
				
			Hashtable<String,Double> targetTable = drugTable.getTargetTable();
			Hashtable<String,Double> pValueTable = drugTable.getProteinNoisePValue();
			
			for(Score score: scores){
				String proteinID = score.proteinID;
				out.write(proteinID+"\t");
					
				//corrDCAVG
				double corrDCAVG = targetTable.get(proteinID);
				if( (corrDCAVG >= -(1+DOUBLE_BUFFER)) && (corrDCAVG <= (1.0+DOUBLE_BUFFER)) ){
					out.write(corrDCAVG+"\t");
				}else{
					out.write("N/A\t");
				}
					
				//PValue and FDR
				double pValue = pValueTable.get(proteinID);
				double FDRvalue = score.FDR;
				
				if( (pValue > (0-DOUBLE_BUFFER)) && (pValue <= (1.0+DOUBLE_BUFFER))){
					if( (FDRvalue > (0-DOUBLE_BUFFER)) && (FDRvalue <= (1.0+DOUBLE_BUFFER))){
						out.write(pValue+"\t"+FDRvalue+"\t");		
					} else {
						out.write(pValue+"\t"+"N/A\t");
					}
				} else{
					for(int j = 0; j < 2; j++){
						out.write("N/A\t");
					}
				}
					
				//prints fold change value
				if( ( foldChange.get(proteinID) >= (0-DOUBLE_BUFFER)) ){
					out.write(foldChange.get(proteinID)+"\t");
				} else{
					out.write("N/A\t");
				}

				// prints if it is significant or not
				Boolean isSignificant = (foldChange.get(proteinID) >= fold_change_threshold) && (FDRvalue > -(0+DOUBLE_BUFFER)) && (FDRvalue <= FDR_threshold);
				out.write(isSignificant+"\n");
				out.flush();
			}
			out.write("\n");
			out.flush();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	
	}

}
