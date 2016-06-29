import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Hashtable;


public class ResultTable {
	// to resolve sticky bit double errors
	private static final double DOUBLE_BUFFER = 0.00000001;
	
	/*
	 * Stores the number of times a certain distribution comes up. Currently
	 * size of 201. 0 corresponds to -1.00 100 corresponds to 0.00 201
	 * corresponds to 1.00
	 */
	private double[][] frequencyMatrix; 
	private double[][] smoothedFrequencyMatrix;
	private double[][] probabilityMatrix;
	private double[] averageDistribution;
	private double[] averageProbability;
	private Hashtable<String, double[]> proteinNoiseProbability;
	private Hashtable<String, Double> proteinNoisePValue;
	private Hashtable<String, double[]> excludeTable;
	private Hashtable<String, double[]> includeTable;
	private Hashtable<String, Double> targetTable;
	private double[][] falsePositiveCounter;
	private int numValidPValues;
	private boolean isDrugTable;
	private static Hashtable<String, Double> foldChange;
	private static Hashtable<String, String> proteinIDtoDescription;
	
	private double foldChangeThreshold;
	
	public ResultTable(boolean _isDrugTable){
		isDrugTable = _isDrugTable;
		numValidPValues = 0;
	}
	
	public void setFrequencyMatrix (double[][] _frequencyMatrix){
		frequencyMatrix = _frequencyMatrix ;
	}
	
	public void setSmoothedFrequencyMatrix (double[][] _smoothedFrequencyMatrix){
		smoothedFrequencyMatrix = _smoothedFrequencyMatrix  ;
	}

	public void setProbabilityMatrix (double[][] _probabilityMatrix){
		probabilityMatrix = _probabilityMatrix;
	}

	public void setAverageDistribution (double[] _averageDistribution){
		averageDistribution = _averageDistribution;
	}

	public void setAverageProbability (double[] _averageProbability){
		averageProbability = _averageProbability;
	}
	
	public void  setProteinNoiseProbability (Hashtable<String, double[]> _proteinNoiseProbability){
		proteinNoiseProbability = _proteinNoiseProbability;
	}

	public void  setProteinNoisePValue (Hashtable<String, Double> _proteinNoisePValue){
		proteinNoisePValue = _proteinNoisePValue;
		countValidPValue();
	}

	public void  setExcludeTable (Hashtable<String, double[]> _excludeTable){
		excludeTable = _excludeTable;
	}

	public void  setIncludeTable (Hashtable<String, double[]> _includeTable){
		includeTable = _includeTable;
	}

	public void  setTargetTable (Hashtable<String, Double> _targetTable){
		targetTable = _targetTable;
	}

	public void setFalsePositiveCounter(double[][] _falsePositiveCounter) {
		falsePositiveCounter = _falsePositiveCounter;
	}

	public void setFoldChange(Hashtable<String, Double> _foldChange, double _foldChangeThreshold){
		foldChange = _foldChange;
		foldChangeThreshold = _foldChangeThreshold;
	}

	public void setProteinIDtoDescription(Hashtable<String, String> _proteinIDtoDescription){
		proteinIDtoDescription = _proteinIDtoDescription;
	}
	
	public double[][] getFrequencyMatrix (){
		return frequencyMatrix;
	}
	
	public double[][] getSmoothedFrequencyMatrix (){
		return smoothedFrequencyMatrix;
	}

	public double[][] getProbabilityMatrix (){
		return probabilityMatrix;
	}

	public double[] getAverageDistribution (){
		return averageDistribution;
	}

	public double[] getAverageProbability (){
		return averageProbability;
	}
	
	public Hashtable<String, double[]> getProteinNoiseProbability (){
		return proteinNoiseProbability;
	}

	public Hashtable<String, Double> getProteinNoisePValue (){
		return proteinNoisePValue;
	}
	
	public Hashtable<String, double[]> getExcludeTable(){
		return excludeTable;
	}
	
	public Hashtable<String, double[]> getIncludeTable(){
		return includeTable;
	}

	public Hashtable<String, Double> getTargetTable(){
		return targetTable;
	}
	
	public double[][] getFalsePositiveCounter() {
		return falsePositiveCounter;
	}
	
	public int getNumValidPValues(){
		return numValidPValues;
	}
	
	public void printFrequencyMatrix(BufferedWriter out){
		try {
			out.write("Frequency Matrix:\n");
			for (int i = 0; i < frequencyMatrix.length; i++) {
				for (int j = 0; j < frequencyMatrix[i].length; j++) {
					out.write(frequencyMatrix[i][j] + "\t");
				}
				out.write("\n");
				out.flush();
			}
			out.write("\n");
			out.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void printSmoothedFrequencyMatrix(BufferedWriter out){
		try {
			out.write("Smoothed Frequency Matrix:\n");
			for (int i = 0; i < smoothedFrequencyMatrix.length; i++) {
				for (int j = 0; j < smoothedFrequencyMatrix[i].length; j++) {
					if(smoothedFrequencyMatrix[i][j] < .000000000000000000001){
						System.out.println("SOMETHING IS WRONG");
					}
					out.write(smoothedFrequencyMatrix[i][j] + "\t");
				}
				out.write("\n");
				out.flush();
			}
			out.write("\n");
			out.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void printProbabilityMatrix(BufferedWriter out){
		try {
			out.write("Probability Matrix:\n");
			for (int i = 0; i < probabilityMatrix.length; i++) {
				for (int j = 0; j < probabilityMatrix[i].length; j++) {
					out.write(probabilityMatrix[i][j] + "\t");
				}
				out.write("\n");
				out.flush();
			}
			out.write("\n");
			out.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void printAverageDistribution(BufferedWriter out){
		try {
			out.write("Distribution of Average Correlations:\n");
			for (int i = 0; i < averageDistribution.length; i++) {
				out.write(averageDistribution[i] + "\t");
			}
			out.write("\n\n");
			out.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}	
	
	public void printAverageProbability(BufferedWriter out){
		try {
			out.write("Probability of Average Correlations:\n");
			for (int i = 0; i < averageProbability.length; i++) {
				out.write(averageProbability[i] + "\t");
			}
			out.write("\n\n");
			out.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void printProteinNoiseTable(BufferedWriter out){
		//This is for everything after the header
		try {
			for (String proteinID : proteinNoiseProbability.keySet()) {
				String output = "N\t"+proteinID + "\t" + proteinIDtoDescription.get(proteinID) + "\t";
				
				// keeps track of whether to print the fold change or N/A
				boolean printFoldChange = true;
				
				// writes the corr_cont_x_y values to the table
				double sum = 0;
				int count = 0;
				for(double d: includeTable.get(proteinID)){
					if((d >= -(1.0+DOUBLE_BUFFER)) && (d <= (1.0+DOUBLE_BUFFER))){
						output += d+"\t";
						sum += d;
						count++;
					} else{
						output += "N/A\t";
					}
				}
				

				double Corr_cont_avg = (sum/count);
				double Corr_D_avg_C_avg = targetTable.get(proteinID);
				
				if(isDrugTable){
					if( count != 0){
						output += Corr_cont_avg+"\t";
					}else{
						output += "N/A\t";
						//prints N/A for fold change when Corr_cont_avg is N/A
						printFoldChange = false;
					}
				}
				
				// writes the target value to the table
				if(( Corr_D_avg_C_avg >= -(1.0+DOUBLE_BUFFER)) && (Corr_D_avg_C_avg <= (1.0+DOUBLE_BUFFER))){
					output += Corr_D_avg_C_avg+"\t";
				} else{
					output += "N/A\t";
					//prints N/A for fold change when Corr_D_avg_C_avg is N/A
					printFoldChange = false;
				}
				
				if(isDrugTable){
					//TODO
					if( printFoldChange ){
						double value = foldChange.get(proteinID);
						if( value >= (0+DOUBLE_BUFFER) ){
							output += value+"\t";
						} else{
							output += "N/A\t";
						}
					}else{
						output += "N/A\t";
					}
					
					if( printFoldChange ){
						output += (foldChange.get(proteinID) >= foldChangeThreshold)+"\t";
					}else{
						output += "N/A\t";
					}
				}
				
				out.write(output);

				//adds pvalue to the table
				if (proteinNoisePValue.get(proteinID) >= (0-DOUBLE_BUFFER)) {
					out.write(proteinNoisePValue.get(proteinID) + "\t");
				} else {
					out.write("N/A\t");
				}
				
				if(isDrugTable){
					//adds the FPR for the drug table
					double pValue = proteinNoisePValue.get(proteinID);
					boolean hasPrinted = false;
					
					if( (pValue >= -(0+DOUBLE_BUFFER)) && (pValue <= (1.0+DOUBLE_BUFFER))){
						
						if( (pValue > (0-DOUBLE_BUFFER)) && (pValue < DOUBLE_BUFFER) ){
							for(int j = 1; j < falsePositiveCounter[falsePositiveCounter.length-1].length; j++){
								out.write(falsePositiveCounter[falsePositiveCounter.length-1][j]+"\t");
							}
							hasPrinted = true;
						} else{
							// remember FPC goes from 1 down to 0
							for(int i = 0; i < falsePositiveCounter.length-1; i++){
								// this rounds them up. This is okay because you should not have a 0 pValue
								if( (pValue <= falsePositiveCounter[i][0]+DOUBLE_BUFFER) && (pValue > falsePositiveCounter[i+1][0]+DOUBLE_BUFFER) ){
									//prints everything in the FPC but the column with the PValue bin
									for(int j = 1; j < falsePositiveCounter[i].length; j++){
										out.write(falsePositiveCounter[i][j]+"\t");
									}
									hasPrinted = true;
								}
							}	
						}
						
					} else{
						for(int j = 1; j < falsePositiveCounter[0].length; j++){
							out.write("N/A\t");
						}
						hasPrinted = true;
					}
					
					if(!hasPrinted){
						System.err.println(proteinID+" did not print FDR");
					}
					
				}
				
				out.write("*\t");
				

				// adds the noise probability to the table
				for (double d : proteinNoiseProbability.get(proteinID)) {
					if (d >= -(0+DOUBLE_BUFFER)) {
						out.write(d + "\t");
					} else {
						out.write("N/A\t");
					}
				}
				out.write("\n");
			}
			out.write("\n");
		} catch (IOException e){
			e.printStackTrace();
		}

	}
	
	private void countValidPValue(){
		for(String proteinID: proteinNoisePValue.keySet()){
			double d = proteinNoisePValue.get(proteinID);
			if( (d >= -(0+DOUBLE_BUFFER)) && (d <= 1+DOUBLE_BUFFER) ){
				numValidPValues++;
			}
		}
	}
	
	
}
