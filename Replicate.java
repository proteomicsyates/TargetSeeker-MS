import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;

public class Replicate {
	private static final String CONTROL_CONDITION = "C";
	private static final String DRUG_CONDITION = "D";
	private static final String REPLICATE_TAG = "R";
	private static final String FRACTION_TAG = "F";
	// to resolve sticky bit double errors
	private static final double DOUBLE_BUFFER = 0.00000001;

	private static int maxNumberOfSpecCounts;
	private static int numFracs;
	private static int specCountThreshold;
	private String conditionNumber;
	private int replicateNumber;
	private Hashtable<String, List<Double>> ProteinIDtoSpecCountVectors;
	private Hashtable<String, List<Double>> ProteinIDtoSpecCountVectorsNormalizedAcrossFractions;
	private Hashtable<String, List<Double>> ProteinIDtoSpecCountVectorsNormalizedAcrossProteinID;
	// the order in which we add the spec count (sorts it)
	private int[] fracOrder;
	// counts how many spectral counts total there are for each rep
	private Hashtable<String, Double> specAbundance;
	// counts how many fracs in each rep are higher than 0 spec counts
	private Hashtable<String, Integer> fracAbundance;

	public Replicate(String _conditionNumber, int _replicateNumber) {
		conditionNumber = _conditionNumber;
		replicateNumber = _replicateNumber;
		fracOrder = null;

		ProteinIDtoSpecCountVectors = new Hashtable<String, List<Double>>();
	}

	// sets the spec count threshold
	public void setSpecCountThreshold(int thresh) {
		specCountThreshold = thresh;
		// expected to be a whole number > 0
		// That is, spec counts that are less than this threshold will not be
		// normalized
		// and will be left out of analysis.
	}

	// maps a spec count vector to a protein ID
	public void loadSpecCounts(String proteinID, List<Double> SpecCountVector) {
		ProteinIDtoSpecCountVectors.put(proteinID, SpecCountVector);
	}

	// Finds and stores the order in which to search the file for each rep
	// Searches the header of the file for the exact column name in the
	// following format:
	// CRxFy for control and DRxFy for drug
	public Boolean setOrder(String[] specHeader, int _numFracs) {
		numFracs = _numFracs;
		fracOrder = new int[numFracs];

		// Looks through the header for the position of each fraction
		for (int i = 0; i < numFracs; i++) {
			Boolean columnFound = false;
			// the column name we are searching for
			String workingName = toString() + FRACTION_TAG + (i + 1);
			// searches through each position in the header
			for (int j = 1; j < specHeader.length; j++) {
				String workingHeaderName = specHeader[j];
				if (workingName.equals(workingHeaderName)) {
					fracOrder[i] = j;
					columnFound = true;
					break;
				}
			}
			if (!columnFound) {
				return columnFound;
			}
		}
		return true;
	}

	// takes in a frac order
	public void setOrder(int[] order) {
		fracOrder = order;
	}

	public int getRep() {
		return replicateNumber;
	}

	public String getCondition() {
		return "Condition " + conditionNumber;
	}

	// returns the indeces in which to gather its fractions
	public int[] fracOrder() {
		return fracOrder;
	}

	// returns the spec count, separated by tabs for easy printing
	public String getSpecForOutput(String proteinID) {
		List<Double> SpecCountVector = ProteinIDtoSpecCountVectors.get(proteinID);
		String output = new String();
		if (SpecCountVector != null) {
			for (int i = 0; i < SpecCountVector.size(); i++) {
				output += (SpecCountVector.get(i)) + "\t";
			}
			return output;
		}
		return "FILE NOT FOUND";
	}

	// returns the spec count vector
	public List<Double> getSpecVector(String proteinID) {
		return ProteinIDtoSpecCountVectors.get(proteinID);
	}

	public String getHash() {
		return ProteinIDtoSpecCountVectors.toString();
	}

	public String toString() {
		// of the format CRx for control DRy for drug
		if (conditionNumber.equals(CONTROL_CONDITION)) {
			return CONTROL_CONDITION + REPLICATE_TAG + replicateNumber;
		} else {
			return DRUG_CONDITION + REPLICATE_TAG + replicateNumber;
		}
	}

	public Hashtable<String, List<Double>> getProteinIDtoSpecCountVectors() {
		return ProteinIDtoSpecCountVectors;
	}

	public void setAbundance() {
		// holds the max number of total spec counts there are for all reps
		maxNumberOfSpecCounts = 0;
		specAbundance = new Hashtable<String, Double>();
		fracAbundance = new Hashtable<String, Integer>();
		for (String proteinID : ProteinIDtoSpecCountVectors.keySet()) {
			List<Double> specVector = ProteinIDtoSpecCountVectors.get(proteinID);
			double abundance = 0;
			int numFracs = 0;
			for (double d : specVector) {
				abundance += d;
				// excludes 0's in fracAbundance
				if (d > 0) {
					numFracs++;
				}
			}
			specAbundance.put(proteinID, abundance);
			fracAbundance.put(proteinID, numFracs);

			if (abundance > maxNumberOfSpecCounts) {
				maxNumberOfSpecCounts = (int) Math.round(abundance);
			}
		}
	}

	public Hashtable<String, Double> getSpecAbundance() {
		return specAbundance;
	}

	public Hashtable<String, Integer> getFracAbundance() {
		return fracAbundance;
	}

	public static int getMaxNumberOfSpecCounts() {
		return maxNumberOfSpecCounts;
	}

	public void normalizeSpecCountAcrossFractions() {
		ProteinIDtoSpecCountVectorsNormalizedAcrossFractions = new Hashtable<String, List<Double>>();
		Set<String> proteinIDSet = ProteinIDtoSpecCountVectors.keySet();
		double[] sum = new double[numFracs];

		// sums up the spec counts of each fraction
		for (int i = 0; i < numFracs; i++) {
			// for each proteinID, adds the spec count at index i
			for (String proteinID : proteinIDSet) {
				sum[i] += ProteinIDtoSpecCountVectors.get(proteinID).get(i);
			}
		}

		// normalizes and adds the spec vectors for each proteinID
		for (String proteinID : proteinIDSet) {
			List<Double> SpecCountVector = ProteinIDtoSpecCountVectors.get(proteinID);
			List<Double> NormalizedSpecCountVector = new ArrayList<Double>();

			for (int i = 0; i < numFracs; i++) {
				if (sum[i] > (0 + DOUBLE_BUFFER)) { // excludes 0's. sum
													// should be a whole
					// number
					NormalizedSpecCountVector.add(SpecCountVector.get(i) / sum[i]);
				} else {
					// probably will not happen, but good to be safe
					NormalizedSpecCountVector.add(0.0);
				}
			}
			ProteinIDtoSpecCountVectorsNormalizedAcrossFractions.put(proteinID, NormalizedSpecCountVector);
		}

	}

	public void normalizeSpecCountAcrossProteinID() {
		ProteinIDtoSpecCountVectorsNormalizedAcrossProteinID = new Hashtable<String, List<Double>>();
		Boolean valuesLessThanOne = true;

		// check if values are less than one -> don't apply spec count threshold
		lessThanOneLoop: for (int i = 0; i < numFracs; i++) {
			for (String proteinID : ProteinIDtoSpecCountVectors.keySet()) {
				if (ProteinIDtoSpecCountVectors.get(proteinID).get(i) > (1 + DOUBLE_BUFFER)) {
					valuesLessThanOne = false;
					break lessThanOneLoop;
				}
			}
		}
		for (String proteinID : ProteinIDtoSpecCountVectors.keySet()) {
			List<Double> SpecCountVector = ProteinIDtoSpecCountVectors.get(proteinID);
			List<Double> NormalizedSpecCountVector = new ArrayList<Double>();

			// sum up all original spec counts
			double sum = 0;
			for (int i = 0; i < SpecCountVector.size(); i++) {
				sum += SpecCountVector.get(i);
			}

			// normalizes the spec counts across the proteinID
			if (((sum > (specCountThreshold + DOUBLE_BUFFER)) && (!valuesLessThanOne))
					|| ((sum > (0 + DOUBLE_BUFFER) && (valuesLessThanOne)))) {
				// bigger than 5.
				// so 6 and up.
				for (int i = 0; i < SpecCountVector.size(); i++) {
					NormalizedSpecCountVector.add(SpecCountVector.get(i) / sum);
				}
			} else {
				for (int i = 0; i < SpecCountVector.size(); i++) {
					NormalizedSpecCountVector.add(-2.0);
				}
			}

			ProteinIDtoSpecCountVectorsNormalizedAcrossProteinID.put(proteinID, NormalizedSpecCountVector);
		}

		// for (String proteinID :
		// ProteinIDtoSpecCountVectorsNormalizedAcrossProteinID.keySet()) {
		// List<Double> SpecCountVector =
		// ProteinIDtoSpecCountVectorsNormalizedAcrossProteinID.get(proteinID);
		// String output = proteinID;
		// for (Double d : SpecCountVector) {
		// output += "\t" + d;
		// }
		// System.out.println(output);
		// }
	}

	public Hashtable<String, List<Double>> getProteinIDtoSpecCountVectorsNormalizedAcrossFractions() {
		return ProteinIDtoSpecCountVectorsNormalizedAcrossFractions;
	}

	public Hashtable<String, List<Double>> getProteinIDtoSpecCountVectorsNormalizedAcrossProteinID() {
		return ProteinIDtoSpecCountVectorsNormalizedAcrossProteinID;
	}
}
