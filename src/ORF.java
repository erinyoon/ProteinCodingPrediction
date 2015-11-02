public class ORF {
	int start;
	int finish;
	int length;
	boolean isCDS = false;
	boolean markov_isCDS = false;
	double markov_score = 0;
	
	public ORF(int start, int finish) {
		this.start = start;
		this.finish = finish;
		this.length = (finish - start + 1);
	}

	public String toString() {
		return start + "\t" + length + "\t" + isCDS + "\t" + markov_score;
	}
}