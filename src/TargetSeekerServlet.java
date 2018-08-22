import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Random;

import javax.servlet.RequestDispatcher;
import javax.servlet.ServletException;
import javax.servlet.annotation.MultipartConfig;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

@javax.servlet.annotation.MultipartConfig
public class TargetSeekerServlet extends HttpServlet {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String filePath;

	@Override
	public void init() {
		// Get the file location where it would be stored.
		// filePath =
		// getServletContext().getInitParameter("file-upload");

	}

	@Override
	public void doPost(HttpServletRequest request, HttpServletResponse response)
			throws ServletException, java.io.IOException {

		System.out.println("Hello!!!!");
		BufferedWriter logD = null;

		try {
			logD = new BufferedWriter(
					new FileWriter("/home/target/tomcat/webapps/Target_Seeker/output/LogDebug.txt", true));
			logD.write("VERY Start DoPost");
			logD.newLine();
			logD.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}

		try {

			boolean isMultipart;
			File file;

			int numControlReplicates = 0;
			int numDrugReplicates = 0;
			int numOfFractions = 0;
			String filePathOnServer = null;
			int kValue = 20;
			double FDRThreshold = 0.1;
			double foldChangeThreshold = 0.2;
			String email = null;

			// Check that we have a file upload request
			java.io.PrintWriter out = response.getWriter();

			try {
				logD.write("afterGetWriter");
				logD.newLine();
				logD.flush();
			} catch (Exception e) {
				e.printStackTrace();
			}

			numControlReplicates = Integer.parseInt(request.getParameter("ControlReplicates"));
			numDrugReplicates = Integer.parseInt(request.getParameter("DrugReplicates"));
			numOfFractions = Integer.parseInt(request.getParameter("fractions"));
			kValue = Integer.parseInt(request.getParameter("kValue"));
			FDRThreshold = Double.parseDouble(request.getParameter("FDRThreshold"));
			foldChangeThreshold = Double.parseDouble(request.getParameter("FoldChange"));
			email = request.getParameter("email");

			try {
				logD.write("afterGettingParams");
				logD.newLine();
				logD.flush();
			} catch (Exception e) {
				e.printStackTrace();
			}

			javax.servlet.http.Part filePart = request.getPart("file");
			String filename = getFilename(filePart);
			InputStream filecontent = filePart.getInputStream();
			OutputStream outputStream = null;

			try {
				logD.write("afterRequestingFile");
				logD.newLine();
				logD.flush();
			} catch (Exception e) {
				e.printStackTrace();
			}

			Random r = new Random();
			String newFileName = filename + r.nextInt();
			filePathOnServer = "/home/target/tomcat/webapps/data/" + newFileName;

			int lengthOfFileExtension = 9;
			String leadingZeros = "";
			for (int i = 0; i < lengthOfFileExtension; i++) {
				leadingZeros += "0";
			}

			File f = new File(filePathOnServer);
			while (f.exists() && !f.isDirectory()) {
				Integer i = (int) (r.nextDouble() * Math.pow(10, lengthOfFileExtension));
				String fileExtension = leadingZeros.substring(0, lengthOfFileExtension - (i.toString()).length())
						+ i.toString();
				newFileName = filename + fileExtension;
				filePathOnServer = "/home/mlaval/tomcat/webapps/data/" + newFileName;
				f = new File(filePathOnServer);
			}

			try {
				logD.write("BeforeLog");
				logD.newLine();
				logD.flush();
			} catch (Exception e) {
				e.printStackTrace();
			} finally {
				if (logD != null)
					try {
						logD.close();
					} catch (Exception e) {
						e.printStackTrace();
					}
			}

			outputStream = new FileOutputStream(new File(filePathOnServer));

			// Log for uploaded files
			DateFormat df = new SimpleDateFormat("dd/MM/yy HH:mm:ss");
			Date dateobj = new Date();

			BufferedWriter log = null;

			try {
				log = new BufferedWriter(
						new FileWriter("/home/target/tomcat/webapps/Target_Seeker/output/LogUsage.txt", true));
				log.write(email + "\t" + newFileName + "\t" + df.format(dateobj) + "\tnumControlReplicates:"
						+ numControlReplicates + "\tnumDrugReplicates:" + numDrugReplicates + "\tnumOfFractions:"
						+ numOfFractions + "\tkValue:" + kValue + "\tFDRThreshold:" + FDRThreshold
						+ "\tfoldChangeThreshold:" + foldChangeThreshold);
				log.newLine();
				log.flush();
			} catch (Exception e) {
				e.printStackTrace();
			} finally {
				if (log != null)
					try {
						log.close();
					} catch (Exception e) {
						e.printStackTrace();
					}
			}

			int read = 0;
			byte[] bytes = new byte[1024];

			while ((read = filecontent.read(bytes)) != -1) {
				outputStream.write(bytes, 0, read);
			}
			System.out.println("Done!");
			if (filecontent != null) {
				filecontent.close();
				if (outputStream != null) {
					outputStream.close();
				}
			}

			// isMultipart = ServletFileUpload.isMultipartContent(request);
			response.setContentType("text/html");
			BufferedWriter outNewFile = new BufferedWriter(new FileWriter(
					new File("/home/target/tomcat/webapps/Target_Seeker/output/" + newFileName + ".txt")));
			outNewFile.write(
					"Your Target-Seeker-MS results are being computed. This may take a few minutes. Please wait...");
			outNewFile.flush();
			// Put a GIF of a spinning wheel for waiting.
			outNewFile.write("The results will be available shortly at this URL");
			outNewFile.flush();
			outNewFile.close();
			String downloadAddress = "http://sealion.scripps.edu:28080/Target_Seeker/output/" + newFileName + ".txt";

			/*
			 * out.println("<html>"); out.println("<head>"); out.println(
			 * "<title>Servlet upload</title>"); out.println("</head>");
			 * out.println("<body>");
			 * 
			 * out.println(
			 * "Your Target-Seeker-MS results are being computed. This may take a few minutes. Please wait..."
			 * ); out.println("</br>"); out.println("</br>");
			 * 
			 * out.println("Parameters: </br>"); out.println(
			 * "Number of control replicates: " + numControlReplicates +
			 * "</br>"); out.println("Number of drug replicates: " +
			 * numDrugReplicates + "</br>"); out.println("Number of fractions: "
			 * + numOfFractions+ "</br>"); out.println("k-value: " + kValue+
			 * "</br>"); out.println("FDR threshold: " + FDRThreshold+ "</br>");
			 * out.println("Fold-change threshold: " + foldChangeThreshold+
			 * "</br>"); out.println("</br>"); out.println("</br>");
			 * out.println(
			 * "<img src=\"http://sealion.scripps.edu:28080/Target_Seeker/loading.GIF\">"
			 * ); out.println("</br>"); out.println("</br>"); out.println(
			 * "You will be able to download your result file at this address: <a href="
			 * + downloadAddress +">"+downloadAddress+"</a>");
			 * out.println("</br>"); out.println("</body>");
			 * out.println("</html>");
			 * 
			 * out.flush(); out.close();
			 */

			String outputFileName = newFileName + ".txt";

			String outputPathAndFile = "/home/target/tomcat/webapps/Target_Seeker/output/" + newFileName + ".txt";
			// Create object of Main class
			// Run the beginning of the program (formally, the main function)
			TargetSeekerMS targetSeekerMS = new TargetSeekerMS(filePathOnServer, numControlReplicates,
					numDrugReplicates, numOfFractions, kValue, FDRThreshold, foldChangeThreshold, outputPathAndFile);
			targetSeekerMS.run();

			/*
			 * try { Thread.sleep(10000); //1000 milliseconds is one second. }
			 * catch(InterruptedException ex) {
			 * Thread.currentThread().interrupt(); }
			 */

			String nextJSP = "/output/" + newFileName + ".txt";
			RequestDispatcher dispatcher = getServletContext().getRequestDispatcher(nextJSP);
			dispatcher.forward(request, response);

			System.out.println("REALLY_DONE");
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	@Override
	public void doGet(HttpServletRequest request, HttpServletResponse response)
			throws ServletException, java.io.IOException {

	}

	public static String getFilename(javax.servlet.http.Part part) {
		for (String cd : part.getHeader("content-disposition").split(";")) {
			if (cd.trim().startsWith("filename")) {
				String filename = cd.substring(cd.indexOf('=') + 1).trim().replace("\"", "");
				return filename.substring(filename.lastIndexOf('/') + 1).substring(filename.lastIndexOf('\\') + 1);
			}
		}
		return null;
	}

}
