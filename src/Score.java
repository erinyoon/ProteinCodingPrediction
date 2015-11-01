
public class Score {
	int A;
	int C;
	int G;
	int T;
	int sum;
	
	public Score(int a, int c, int g, int t, int count) {
		A = a; C = c; G = g; T = t; sum = count;
	}
	
	public int getTotal() {
		return A + C + G + T;
	}
	public String toString() {
		return "A:" + A + " C:" + C + " G:" + G + " T:" + T;
	}
}
