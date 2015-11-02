
public class Score {
	int A; // each seen
	int C;
	int G;
	int T;
	
	public Score(int a, int c, int g, int t) {
		A = a; C = c; G = g; T = t;
	}
	
	public int getTotal() {
		return A + C + G + T;
	}
	
	public double getProb(String temp) {
		if (temp.equalsIgnoreCase("a")) {
			return 1.0 * A / getTotal();
		} else if (temp.equalsIgnoreCase("c")) {
			return 1.0 * C / getTotal();
		} else if (temp.equalsIgnoreCase("g")) {
			return 1.0 * G / getTotal();
		} else {
			return 1.0 * T / getTotal();
		}
	}
	
	public String toString() {
		return "A:" + A + " C:" + C + " G:" + G + " T:" + T;
	}
}
