import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

public class TargetSeekerMSCommandLine {
	private static Options options;

	/*** Default parameters ***/
	public static final double DEFAULT_PPM_TOLERANCE = 0.000002;
	public static final double DEFAULT_X_CORR_SCORE = 2.0;
	public static final int DEFAULT_MISSED_CLEAVAGES = 0;
	public static final ExclusionProfileEnum DEFAULT_EXCLUSION_PROFILE = ExclusionProfileEnum.MACHINE_LEARNING_GUIDED_EXCLUSION_PROFILE;

	/*** File paths ***/
	private static String fasta_file_name;
	private static String digested_fasta_file_name;
	private static String mzml_file_name;
	private static String mzid_file_name;
	private static String result_database_file_name;

	/*** Analysis parameters ***/
	private static double ppmTolerance;
	private static int missedCleavages;
	// threshold for xCorr confidence
	private static double xCorrScore;
	// If this protein has confidently been IDd more than 3 times, add it to the
	// protein, observed mass, and theoretical peptide masses to the exclusion list
	private static final int numDB = 3;
	private static ExclusionProfileEnum exclusionProfileType;

	private static Simulator simulator;
	private static ExclusionProfile exclusionProfile;

	public static void main(String[] args) {
		log.info("Starting program");

		// Parse the command line options
		parseArguments(args);
		// Load files necessary for program
		loadFiles();
		// Run simulation
		runSimulation();

		log.debug("End of program.");

	}

	private static void runSimulation() {
		log.debug("Running Real Time Simulation");
		simulator.runSimulation(exclusionProfile);
	}

	private static void loadFiles() {
		Database database;
		ResultDatabase resultDatabase;
		MZMLFile mzml;

		// Load the database file
		if (digested_fasta_file_name == null) {
			log.debug("Loading fasta and digesting the fasta file");
			database = Loader.loadExclusionDatabase(fasta_file_name, missedCleavages);
		} else {
			log.debug("Loading both fasta and digested fasta file");
			database = Loader.loadExclusionDatabase(fasta_file_name, digested_fasta_file_name);
		}

		// Load the mzml file needed for the result database and the simulator
		mzml = Loader.parseMZML(mzml_file_name);

		// Load the result database file
		if (mzid_file_name == null && result_database_file_name == null) {
			log.debug("No results database. Stepping through the simulation blind.");
			resultDatabase = null;
		} else {
			if (result_database_file_name != null) {
				log.debug("Loading the result database");
				resultDatabase = Loader.parseResultDatabase(result_database_file_name);
			} else {
				log.debug("Creating the result database");
				resultDatabase = Loader.loadResultDatabase(mzml, mzid_file_name);
			}
		}

		// Load the simulator
		log.debug("Initializing simulator");
		simulator = new Simulator(mzml.getSpectraArray());

		// Load the Exclusion Profile
		log.debug("Initializing exclusion profile...");
		log.debug("Exclusion profile: " + exclusionProfileType);
		switch (exclusionProfileType) {
		case MACHINE_LEARNING_GUIDED_EXCLUSION_PROFILE:
			exclusionProfile = new MachineLearningGuidedExclusion(null, database); // TODO not done
			break;
		case NORA_EXCLUSION_PROFILE:
			// exclusionProfile = new NoraExclusion(database, resultDatabase, xCorrScore,
			// ppmTolerance, numDB); // TODO not
			// done
			break;
		}

		// simulator.test(database, resultDatabase);
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
			if (cmd.hasOption("d")) {
				// will display debug messages to console
				log.setLevel(Level.DEBUG);
				log.debug("Debug mode activated");
			} else {
				// will not display debug messages to console
				log.setLevel(Level.INFO);
			}

			/*** File paths ***/
			fasta_file_name = cmd.getOptionValue("f");
			digested_fasta_file_name = cmd.getOptionValue("g", null);
			mzml_file_name = cmd.getOptionValue("z");
			mzid_file_name = cmd.getOptionValue("y", null);
			result_database_file_name = cmd.getOptionValue("r", null);
			// checking if file paths are valid
			if (!Loader.checkValidFilePath(fasta_file_name)) {
				throw (new ParseException("Invalid fasta file"));
			} else if ((digested_fasta_file_name != null) && !Loader.checkValidFilePath(digested_fasta_file_name)) {
				throw (new ParseException("Invalid digested fasta file"));
			} else if (!Loader.checkValidFilePath(mzml_file_name)) {
				throw (new ParseException("Invalid mzML file"));
			} else if ((digested_fasta_file_name != null) && !Loader.checkValidFilePath(mzid_file_name)) {
				throw (new ParseException("Invalid mzid file"));
			} else if ((result_database_file_name != null) && !Loader.checkValidFilePath(result_database_file_name)) {
				throw (new ParseException("Invalid ResultDatabase file"));
			}

			/*** Program parameters ***/
			// missed cleavages
			if (cmd.hasOption("n")) { // TODO is the missed cleavages necessary if you include digested file?
				try {
					missedCleavages = Integer.parseInt(cmd.getOptionValue("n"));
				} catch (NumberFormatException e) {
					throw new ParseException("-n must be an integer");
				}
			} else {
				missedCleavages = DEFAULT_MISSED_CLEAVAGES;
			}

			// exclusion profile
			if (cmd.hasOption("e")) {
				try {
					switch (Integer.parseInt(cmd.getOptionValue("e"))) {
					case 0:
						exclusionProfileType = ExclusionProfileEnum.NORA_EXCLUSION_PROFILE;
						break;
					case 1:
						exclusionProfileType = ExclusionProfileEnum.MACHINE_LEARNING_GUIDED_EXCLUSION_PROFILE;
						break;
					default:
						exclusionProfileType = DEFAULT_EXCLUSION_PROFILE;
						break;
					}
				} catch (NumberFormatException e) {
					throw new ParseException("-e must be an integer");
				}
			} else {
				exclusionProfileType = DEFAULT_EXCLUSION_PROFILE;
			}

			// xCorr threshold
			if (cmd.hasOption("x")) {
				try {
					xCorrScore = Double.parseDouble(cmd.getOptionValue("x"));
				} catch (NumberFormatException e) {
					throw new ParseException("-x must be an double");
				}
			} else {
				xCorrScore = DEFAULT_X_CORR_SCORE;
			}

			// ppmTolerance threshold
			if (cmd.hasOption("p")) { // TODO make sure the ppm is a small number, and not like '5' or something.
				try {
					ppmTolerance = Double.parseDouble(cmd.getOptionValue("p"));
				} catch (NumberFormatException e) {
					throw new ParseException("-p must be an double");
				}
			} else {
				ppmTolerance = DEFAULT_PPM_TOLERANCE;
			}

		} catch (ParseException e) {
			errorInArguments(e.getMessage());
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println(e.getMessage());
		}
		log.debug("Fasta file name: " + fasta_file_name);
		log.debug("Digested fasta file name: " + digested_fasta_file_name);
		log.debug("mzML file name: " + mzml_file_name);
		log.debug("mzid file name: " + mzid_file_name);
		log.debug("Result database file name: " + result_database_file_name);
		log.debug("ppmTolerance: " + ppmTolerance);
		log.debug("missedCleavages: " + missedCleavages);
		log.debug("xCorr threshold: " + xCorrScore);
		log.debug("ExclusionProfile: " + exclusionProfileType);
	}

	private static void errorInArguments(String header) {
		// automatically generate the help statement
		HelpFormatter formatter = new HelpFormatter();
		formatter.setOptionComparator(null);
		if (header == null) {
			formatter.printHelp("java -jar RealTimeMS.jar", options);
		} else {
			formatter.printHelp("java -jar RealTimeMS.jar", "\n************\n" + header + "\n************\n", options,
					"Contact Alexander Pelletier at apell035@uottawa.ca for more help");
		}
		System.exit(0);
	}

	private static void setupOptions() {
		// TODO Edit descriptions. Also check if some need to be capitalized, for the
		// options.
		// create Options object
		options = new Options();

		// add option
		options.addOption("d", "debug", false, "Debug mode.");
		options.addOption("e", "exclusionProfile", true, "0 - Nora Exclusion, 1 - MachineLearningGuidedMS (default).");
		options.addOption("f", "fasta", true, "Fasta file to be used as a database.");
		options.addOption("g", "digestedFasta", true, "Digested fasta file to be used as a database.");
		options.addOption("h", "help", false, "Displays help options.");
		options.addOption("n", "missedCleavages", true, "Number of missed cleavages allowed for tryptic digestion.");
		options.addOption("o", "xCorrScore", true, "Cross correlation cuttoff threshold (depreciated)");
		options.addOption("p", "ppmTolerance", true, "ppmTolerance allowed by the analysis.");
		options.addOption("r", "resultDatabase", true, "Result database file (.tsv).");
		options.addOption("y", "mzid", true, "Identification file");
		options.addOption("z", "mzML", true, "MS experiment to be used");

	}

	}

	public static void main(String[] args) {
		String outputFileName = newFileName + ".txt";

		String outputPathAndFile = "/home/target/tomcat/webapps/Target_Seeker/output/" + newFileName + ".txt";
		// Create object of Main class
		// Run the beginning of the program (formally, the main function)
		TargetSeekerMS targetSeekerMS = new TargetSeekerMS(filePathOnServer, numControlReplicates,
				numDrugReplicates, numOfFractions, kValue, FDRThreshold, foldChangeThreshold, outputPathAndFile);
		targetSeekerMS.run();
	}
	
}
