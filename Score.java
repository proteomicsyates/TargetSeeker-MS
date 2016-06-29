public class Score implements Comparable<Score>{
	double FDR;
	String proteinID;

	public Score(double FDR, String proteinID){
		this.FDR = FDR;
		this.proteinID = proteinID;
	}

	@Override
	public int compareTo(Score o){
		return this.FDR < o.FDR ? -1: this.FDR > o.FDR ? 1: 0;
	}
}