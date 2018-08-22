import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class TargetSeekerMSCommandLine {
	private static Options options;

	/*** Default parameters ***/
	public static final int DEFAULT_CONTROL_REPLICATES = 4;
	public static final int DEFAULT_DRUG_REPLICATES = 3;
	public static final int DEFAULT_NUM_FRACTIONS = 10;
	public static final int DEFAULT_KVALUE = 20;
	public static final double DEFAULT_FDR_THRESHOLD = 0.05;
	public static final double DEFAULT_FOLDCHANGE_THRESHOLD = 0.2;
	public static final String DEFAULT_OUTPUT_FILE_NAME = "TargetSeekerResults.txt";

	/*** Analysis parameters ***/
	private static String inputFileName;
	private static int numControlReplicates;
	private static int numDrugReplicates;
	private static int numOfFractions;
	private static int kValue;
	private static double FDRThreshold;
	private static double foldChangeThreshold;
	private static String outputFileName;
	
	public static void main(String[] args) {
		System.out.println("Starting program");
		// Parse the command line options
		parseArguments(args);
		// Run simulation
		execute();
		System.out.println("End of program.");
	}

	private static void execute() {
		System.out.println("Running TargetSeeker-MS");
		// Run the beginning of the program (formally, the main function)
		TargetSeekerMS targetSeekerMS = new TargetSeekerMS(inputFileName, numControlReplicates,
				numDrugReplicates, numOfFractions, kValue, FDRThreshold, foldChangeThreshold, outputFileName);
		targetSeekerMS.run();
	}

	private static void parseArguments(String[] args) {
		// parse args
		setupOptions();
		CommandLineParser parser = new DefaultParser();
		try {
			CommandLine cmd = parser.parse(options, args);

			/*** Help and debug ***/
			if (cmd.hasOption("h")) {
				// stop parsing and display help message
				throw (new ParseException(null));
			}

			/*** File paths ***/
			inputFileName = cmd.getOptionValue("i");
			outputFileName = cmd.getOptionValue("o");
			// checking if file paths are valid
			if (!checkValidFilePath(inputFileName)) {
				throw (new ParseException("Invalid input file!"));
			} else if ((outputFileName != null) && !checkValidFilePath(outputFileName)) {
				throw (new ParseException("Invalid output file!"));
			} 

			/*** Program parameters ***/
			// numControlReplicates
			if (cmd.hasOption("c")) {
				try {
					numControlReplicates = Integer.parseInt(cmd.getOptionValue("c"));
					if (numControlReplicates < 4) {
						throw new NumberFormatException();
					}
				} catch (NumberFormatException e) {
					throw new ParseException("numControlReplicates must an integer greater than or equal to 4");
				}
			} else {
				numControlReplicates = DEFAULT_CONTROL_REPLICATES;
			}
			// numDrugReplicates
			if (cmd.hasOption("d")) {
				try {
					numDrugReplicates = Integer.parseInt(cmd.getOptionValue("d"));
					if (numDrugReplicates <= 0) {
						throw new NumberFormatException();
					}
				} catch (NumberFormatException e) {
					throw new ParseException("numDrugReplicates must an integer greater than 0");
				}
			} else {
				numDrugReplicates = DEFAULT_DRUG_REPLICATES;
			}
			// numOfFractions
			if (cmd.hasOption("f")) {
				try {
					numOfFractions = Integer.parseInt(cmd.getOptionValue("f"));
					if (numOfFractions <= 2) {
						throw new NumberFormatException();
					}
				} catch (NumberFormatException e) {
					throw new ParseException("numOfFractions must an integer greater than 2");
				}
			} else {
				numOfFractions = DEFAULT_NUM_FRACTIONS;
			}
			// kValue
			if (cmd.hasOption("k")) {
				try {
					kValue = Integer.parseInt(cmd.getOptionValue("k"));
					if (kValue <= 0) {
						throw new NumberFormatException();
					}
				} catch (NumberFormatException e) {
					throw new ParseException("kValue must an integer greater than 0");
				}
			} else {
				kValue = DEFAULT_KVALUE;
			}
			// FDRThreshold
			if (cmd.hasOption("y")) {
				try {
					FDRThreshold = Double.parseDouble(cmd.getOptionValue("y"));
					if (FDRThreshold < 0 || FDRThreshold > 1.0) {
						throw new NumberFormatException();
					}
				} catch (NumberFormatException e) {
					throw new ParseException("FDR threshold must be a number between 0 and 1");
				}
			} else {
				FDRThreshold = DEFAULT_FDR_THRESHOLD;
			}
			// foldChangeThreshold
			if (cmd.hasOption("z")) {
				try {
					foldChangeThreshold = Double.parseDouble(cmd.getOptionValue("z"));
					if (foldChangeThreshold < 0) {
						throw new NumberFormatException();
					}
				} catch (NumberFormatException e) {
					throw new ParseException("Fold-change threshold must be a number greater than 0");
				}
			} else {
				foldChangeThreshold = DEFAULT_FOLDCHANGE_THRESHOLD;
			}
		} catch (ParseException e) {
			errorInArguments(e.getMessage());
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println(e.getMessage());
		}
		System.out.println("inputFileName: " + inputFileName);
		System.out.println("outputFileName: " + outputFileName);
		System.out.println("numControlReplicates: " +numControlReplicates);
		System.out.println("numDrugReplicates: " + numDrugReplicates);
		System.out.println("numOfFractions: " + numOfFractions);
		System.out.println("kValue: " + kValue);
		System.out.println("FDRThreshold: " + FDRThreshold);
		System.out.println("foldChangeThreshold: " + foldChangeThreshold);
	}

	public static boolean checkValidFilePath(String file_path) throws Exception {
		return (new File(file_path).exists());
	}
	
	private static void errorInArguments(String header) {
		// automatically generate the help statement
		HelpFormatter formatter = new HelpFormatter();
		formatter.setOptionComparator(null);
		if (header == null) {
			formatter.printHelp("java -jar TargetSeeker-MS.jar", options);
		} else {
			formatter.printHelp("java -jar TargetSeeker-MS.jar", "\n************\n" + header + "\n************\n", options,
					"Contact Alexander Pelletier at apell035@uottawa.ca for more help");
		}
		System.exit(0);
	}

	private static void setupOptions() {
		// add options	
		options = new Options();
		options.addOption("h", "help", false, "Displays help options.");
		options.addOption("c", "numControlReplicates", true, "number_control_replicates: the number of control replicates (integer >= 4). Default = "+DEFAULT_CONTROL_REPLICATES);
		options.addOption("d", "numDrugReplicates", true, "number_drug_replicates: the number of drug replicates (integer > 0). Default = "+DEFAULT_DRUG_REPLICATES);
		options.addOption("f", "numFractions", true, "the number of fractions in each replicate (integer > 2). Default = "+DEFAULT_NUM_FRACTIONS);
		options.addOption("k", "kValue", true, "the k-nearest-neighbor smoothing factor (integer > 0). Default = "+DEFAULT_KVALUE);
		options.addOption("y", "fdrThreshold", true, "the FDR cut-off for determining significant protein. Any proteins with an FDR higher than this threshold are considered not significant (number between 0 and 1). Default = "+DEFAULT_FDR_THRESHOLD);
		options.addOption("z", "foldChangeThreshold", true, "the fold-change threshold for determining significant protein. Any proteins with a fold-change less than this threshold are considered not significant (non-negative number). Default = "+DEFAULT_FOLDCHANGE_THRESHOLD);
		options.addOption("i", "inputFile", true, ".txt file containing spectral counts separated by tabs. File should follow the following naming convention for header: Protein CR1F1 CR1F2 ... DRxFy");
		options.addOption("o", "outputFile", true, "Output the results to this file. Default = "+DEFAULT_OUTPUT_FILE_NAME);
	}

}
