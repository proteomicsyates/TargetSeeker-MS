import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;

public class Correlation {
	// to resolve sticky bit double errors
	private static final double DOUBLE_BUFFER = 0.00000001;
	
	private static final boolean CORR_DC_REP_X_ON = false;
	private static final boolean CORR_CONT_X_Y_ON = true;
	private static final boolean CORR_CONT_AVG_EXCLUDE_ON = true;
	private static final boolean CORR_D_REP_X_CONT_REP_Y_ON = true;

	private static final boolean PRINT_CORR_DC_REP_X_ON = false;
	private static final boolean PRINT_CORR_CONT_X_Y_ON = true;
	private static final boolean PRINT_CORR_CONT_AVG_EXCLUDE_ON = true;
	private static final boolean PRINT_CORR_D_REP_X_CONT_REP_Y_ON = true;

	private static final int AVG = 100;

	/*
	 * Toggle true/false to compute these values or not Note to self:
	 * Corr_Cont_avg_EXCLUDE relies on Corr_Cont_X_Y Compute frequency matrix
	 * relies on both of the above.
	 * 
	 * This toggling function was originally planned to save computation time
	 * for unneccesary computations. But it seems only Corr_DC_Rep_x is
	 * unnecessary, so far. Therefore, there is a toggle for printing, so if the
	 * flag is false, then the value will not be added to the output file
	 */

	private ArrayList<String> header;
	private Hashtable<String, double[]> proteinCorrelationValue;
	private Hashtable<String, Integer> indexTable;

	private int offset;
	/*
	 * offset accounts for anything before the correlation values. That is,
	 * header.size() gives a certain size, but there is an offset for the
	 * double[], which starts at 0, where it is at index 0+offset.
	 */

	private static int numControlReps[];
	private static int numDrugReps[];

	public Correlation(ArrayList<String> headerColumnsAdded, int[] _numControlReps, int[] _numDrugReps) {
		// headerColumnsAdded should contain PLine and Accession + misc.
		offset = headerColumnsAdded.size();
		int count = 0; // makes sure they are added in order
		indexTable = new Hashtable<String, Integer>(); 
		// maps a name to its index
		proteinCorrelationValue = new Hashtable<String, double[]>(); 
		// maps a name to correlation values
		numControlReps = _numControlReps;
		numDrugReps = _numDrugReps;

		if (CORR_DC_REP_X_ON) {
			for (int i = 0; i < Math.min(numControlReps.length, numDrugReps.length); i++) {
				String temp = "Corr_DC_Rep_" + (numControlReps[i]);
				indexTable.put(temp, count++);
				headerColumnsAdded.add(temp);
			}
		}

		if (CORR_CONT_X_Y_ON) {
			for (int i = 0; i < numControlReps.length; i++) {
				for (int j = (i + 1); j < numControlReps.length; j++) {
					String temp = ("Corr_cont_" + (numControlReps[i]) + "_" + (numControlReps[j]));
					indexTable.put(temp, count++);
					headerColumnsAdded.add(temp);
				}
			}
			String avg = "Corr_cont_avg";
			indexTable.put(avg, count++);
			headerColumnsAdded.add(avg);
		}

		if (CORR_CONT_AVG_EXCLUDE_ON) {
			for (int i = 0; i < numControlReps.length; i++) {
				for (int j = (i + 1); j < numControlReps.length; j++) {
					String temp = ("Corr_cont_avg_EXCLUDE_" + (numControlReps[i]) + "_" + (numControlReps[j]));
					indexTable.put(temp, count++);
					headerColumnsAdded.add(temp);
				}
			}
		}

		if (CORR_D_REP_X_CONT_REP_Y_ON) {
			for (int i = 0; i < numDrugReps.length; i++) {
				for (int j = 0; j < numControlReps.length; j++) {
					String temp = ("Corr_DrugRep_" + (numDrugReps[i]) + "_ContRep_" + (numControlReps[j]));
					indexTable.put(temp, count++);
					headerColumnsAdded.add(temp);
				}
				String avg = "Corr_DrugRep_" + (numDrugReps[i]) + "_ContRep_avg";
				indexTable.put(avg, count++);
				headerColumnsAdded.add(avg);
			}
			String averageAverage = "Corr_DrugRep_avg_ContRep_avg";
			indexTable.put(averageAverage, count++);
			headerColumnsAdded.add(averageAverage);
		}

		header = headerColumnsAdded;
	}

	// initializes each correlation array with a protein ID
	public void initializeSet(Set<String> proteinIDSet) {
		for (String proteinID : proteinIDSet) {
			double[] temp = new double[header.size() - offset];
			proteinCorrelationValue.put(proteinID, temp);
		}
	}

	public boolean isCorrDCRep() {
		return CORR_DC_REP_X_ON;
	}

	public boolean isCorrContXY() {
		return CORR_CONT_X_Y_ON;
	}

	public boolean isCorrDRepCRep() {
		return CORR_D_REP_X_CONT_REP_Y_ON;
	}

	// returns an array list which will be printed in the main header
	public ArrayList<String> getHeader() {
		ArrayList<String> outputHeader = new ArrayList<String>();

		for (int i = 0; i < offset; i++) {
			outputHeader.add(header.get(i));
		}

		if (PRINT_CORR_DC_REP_X_ON) {
			for (int i = 0; i < Math.min(numControlReps.length, numDrugReps.length); i++) {
				String temp = "Corr_DC_Rep_" + (numControlReps[i]);
				outputHeader.add(temp);
			}
		}

		if (PRINT_CORR_CONT_X_Y_ON) {
			for (int i = 0; i < numControlReps.length; i++) {
				for (int j = (i + 1); j < numControlReps.length; j++) {
					String temp = ("Corr_cont_" + (numControlReps[i]) + "_" + (numControlReps[j]));
					outputHeader.add(temp);
				}
			}
			String avg = "Corr_cont_avg";
			outputHeader.add(avg);
		}

		if (PRINT_CORR_CONT_AVG_EXCLUDE_ON) {
			for (int i = 0; i < numControlReps.length; i++) {
				for (int j = (i + 1); j < numControlReps.length; j++) {
					String temp = ("Corr_cont_avg_EXCLUDE_" + (numControlReps[i]) + "_" + (numControlReps[j]));
					outputHeader.add(temp);
				}
			}
		}

		if (PRINT_CORR_D_REP_X_CONT_REP_Y_ON) {
			for (int i = 0; i < numDrugReps.length; i++) {
				for (int j = 0; j < numControlReps.length; j++) {
					String temp = ("Corr_DrugRep_" + (numDrugReps[i]) + "_ContRep_" + (numControlReps[j]));
					outputHeader.add(temp);
				}
				String avg = "Corr_DrugRep_" + (numDrugReps[i]) + "_ContRep_avg";
				outputHeader.add(avg);
			}
			String averageAverage = "Corr_DrugRep_avg_ContRep_avg";
			outputHeader.add(averageAverage);
		}

		return outputHeader;
	}

	// returns the mapping of proteinID to its correlation values
	public Hashtable<String, double[]> getTable() {
		return proteinCorrelationValue;
	}

	public Set<String> keySet() {
		return proteinCorrelationValue.keySet();
	}

	// used for looping through and printing all values
	public double getCorrelationValue(String proteinID, int index) {
		return proteinCorrelationValue.get(proteinID)[index];
	}

	// returns the value of corr_DC_rep_(repNum)
	public double getCorrDCRep(String proteinID, int repNum) {
		double retVal = -2;
		if (CORR_DC_REP_X_ON) {
			String tempName = "Corr_DC_Rep_" + repNum;
			retVal = proteinCorrelationValue.get(proteinID)[indexTable.get(tempName)];
		}

		return retVal;
	}

	// returns the value of corr_cont_x_y
	public double getCorrCont(String proteinID, int repX, int repY) {
		double retVal = -2;
		if (CORR_CONT_X_Y_ON) {
			String tempName;
			if (repX == AVG) {
				tempName = "Corr_cont_avg";
			} else {
				tempName = "Corr_cont_" + repX + "_" + repY;
			}

			retVal = proteinCorrelationValue.get(proteinID)[indexTable.get(tempName)];
		}
		return retVal;
	}

	// returns the value of corr_cont_avg_EXCLUDE_x_y
	public double getCorrContExclude(String proteinID, int repX, int repY) {
		double retVal = -2;
		if (CORR_CONT_AVG_EXCLUDE_ON) {
			String tempName = "Corr_cont_avg_EXCLUDE_" + repX + "_" + repY;
			retVal = proteinCorrelationValue.get(proteinID)[indexTable.get(tempName)];
		}
		return retVal;
	}

	// returns the value of corr_DrugRepX_ContRep_Y
	public double getCorrDCRepXRepY(String proteinID, int drugRepX, int contRepY) {
		double retVal = -2;
		if (CORR_D_REP_X_CONT_REP_Y_ON) {
			String tempName;
			if (drugRepX == AVG) {
				tempName = "Corr_DrugRep_avg_ContRep_avg";
			} else if (contRepY == AVG) {
				tempName = "Corr_DrugRep_" + drugRepX + "_ContRep_avg";
			} else {
				tempName = "Corr_DrugRep_" + drugRepX + "_ContRep_" + contRepY;
			}
			retVal = proteinCorrelationValue.get(proteinID)[indexTable.get(tempName)];
		}
		return retVal;
	}

	// adds the value to its corresponding rep number for corr_DC
	public void addCorrDCRep(String proteinID, double value, int repNum) {
		if (CORR_DC_REP_X_ON) {
			double[] tempArray = proteinCorrelationValue.get(proteinID);
			String tempName = "Corr_DC_Rep_" + repNum;
			tempArray[indexTable.get(tempName)] = value;
			proteinCorrelationValue.put(proteinID, tempArray);
		}
	}

	// adds the values for corr_cont and corr_cont_avg_EXCLUDE in the correct
	// spot
	public void addCorrCont(String proteinID, Hashtable<String, Double> correlationTable) {
		if (CORR_CONT_X_Y_ON) {
			
			double[] tempArray = proteinCorrelationValue.get(proteinID);
			int numElements = 0;
			double sum = 0;
			for (int i = 0; i < numControlReps.length; i++) {
				for (int j = (i + 1); j < numControlReps.length; j++) {
					String tempName = ("Corr_cont_" + (numControlReps[i]) + "_" + (numControlReps[j]));
					double corrVal = correlationTable.get(tempName);
					tempArray[indexTable.get(tempName)] = corrVal;
					
					// only includes in the average if it is valid.
					if( corrVal > -(0+DOUBLE_BUFFER) && corrVal < (1+DOUBLE_BUFFER)){
						sum += correlationTable.get(tempName);
						numElements++;
					}
					
				}
			}
			
			if(numElements > 0){
				tempArray[indexTable.get("Corr_cont_avg")] = sum / numElements;
			}else{
				tempArray[indexTable.get("Corr_cont_avg")] = -2.0;
			}
			
			if (CORR_CONT_AVG_EXCLUDE_ON) {
				for (int i = 0; i < numControlReps.length; i++) {
					for (int j = (i + 1); j < numControlReps.length; j++) {
						String tempName = ("Corr_cont_avg_EXCLUDE_" + (numControlReps[i]) + "_" + (numControlReps[j]));
						double exclude = correlationTable
								.get("Corr_cont_" + (numControlReps[i]) + "_" + (numControlReps[j]));

						if( exclude > -(0+DOUBLE_BUFFER) && exclude < (1+DOUBLE_BUFFER)){
							if(numElements > 1){
								// case 1: works like we expect
								tempArray[indexTable.get(tempName)] = (sum - exclude) / (numElements - 1);
							}else{
								// case 2: there is 0 or only 1 element (which is exclude). prints N/A
								tempArray[indexTable.get(tempName)] = -2.0;
							}
						} else{
							// case 3: exclude is N/A and excluding it does nothing (still avg).
							tempArray[indexTable.get(tempName)] = sum / numElements;
						}
						
					}
				}
			}

			proteinCorrelationValue.put(proteinID, tempArray);
		}
	}

	// adds the values for corr_DrugRepX_ContRepY in the correct spot
	public void addCorrDRepCRep(String proteinID, Hashtable<String, Double> correlationTable) {
		if (CORR_D_REP_X_CONT_REP_Y_ON) {
			double[] tempArray = proteinCorrelationValue.get(proteinID);
			double totalSum = 0; // for calculating the average average
			int totalElements = 0;
			
			for (int i = 0; i < numDrugReps.length; i++) {
				int numElements = 0;
				double sum = 0;
				for (int j = 0; j < numControlReps.length; j++) {
					String tempName = ("Corr_DrugRep_" + (numDrugReps[i]) + "_ContRep_" + (numControlReps[j]));
					double input = correlationTable.get(tempName);
					tempArray[indexTable.get(tempName)] = input;
					// only sums valid correlation values
					if( (input >= -(0+DOUBLE_BUFFER)) && (input <= (1.0+DOUBLE_BUFFER)) ){ 
						sum += correlationTable.get(tempName);
						numElements++;
						totalSum += correlationTable.get(tempName);
						totalElements++;
					}
				}
				String avg = "Corr_DrugRep_" + (numDrugReps[i]) + "_ContRep_avg";
				if(numElements > 0){
					tempArray[indexTable.get(avg)] = sum / numElements;
				} else{
					tempArray[indexTable.get(avg)] = -2.0;
				}
			}

			String averageAverage = "Corr_DrugRep_avg_ContRep_avg";
			if(totalElements!=0){
				tempArray[indexTable.get(averageAverage)] = totalSum / totalElements;
			} else{
				tempArray[indexTable.get(averageAverage)] = -2.0;
			}
			
			proteinCorrelationValue.put(proteinID, tempArray);
		}
	}

	// converts correlation value to its corresponding string, with tab
	public String formatCorrelationValue(double correlation) {
		String output = new String();
		Double correlationObject = new Double(correlation);

		if (correlationObject.isNaN()) { // check if NaN (because one vector
											// might contain the exact same
											// value at every entry.
			output += "N/A\t";
		} else if (correlation != -2.0) { // if not -2.0 then write computed
											// correlation
			output += correlation + "\t";
		} else { // if -2.0 then write N/A
			output += "N/A\t";
		}

		return output;
	}

	// for easy printing of entire class
	public String print(String proteinID) {
		String output = new String();

		if (PRINT_CORR_DC_REP_X_ON) {
			for (int i = 0; i < Math.min(numControlReps.length, numDrugReps.length); i++) {
				output += formatCorrelationValue(getCorrDCRep(proteinID, (numControlReps[i])));
			}
		}

		if (PRINT_CORR_CONT_X_Y_ON) {
			for (int i = 0; i < numControlReps.length; i++) {
				for (int j = (i + 1); j < numControlReps.length; j++) {
					output += formatCorrelationValue(getCorrCont(proteinID, (numControlReps[i]), (numControlReps[j])));
				}
			}
			output += formatCorrelationValue(getCorrCont(proteinID, AVG, AVG));
		}

		if (PRINT_CORR_CONT_AVG_EXCLUDE_ON) {
			for (int i = 0; i < numControlReps.length; i++) {
				for (int j = (i + 1); j < numControlReps.length; j++) {
					output += formatCorrelationValue(
							getCorrContExclude(proteinID, (numControlReps[i]), (numControlReps[j])));
				}
			}
		}

		if (PRINT_CORR_D_REP_X_CONT_REP_Y_ON) {
			for (int i = 0; i < numDrugReps.length; i++) {
				for (int j = 0; j < numControlReps.length; j++) {
					output += formatCorrelationValue(
							getCorrDCRepXRepY(proteinID, (numDrugReps[i]), (numControlReps[j])));
				}
				output += formatCorrelationValue(getCorrDCRepXRepY(proteinID, (numDrugReps[i]), AVG));
			}
			output += formatCorrelationValue(getCorrDCRepXRepY(proteinID, AVG, AVG));
		}

		return output;
	}
}
