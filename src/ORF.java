
public class ORF {
	int start;
	int finish;
	int count;
	boolean isCDS = false;
	double markov = 0.0;
	
	public ORF(int start, int finish) {
		this.start = start;
		this.finish = finish;
		this.count = (finish - start + 1);
	}
	
	public String toString() {
		return "[" + start + "," + finish + "]";
	}
}
