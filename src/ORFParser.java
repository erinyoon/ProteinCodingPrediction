import java.util.HashSet;


/**
 * ORFParser takes in file of Methanocaldococcus jannaschii
 *  
 *  NOTE:
 *  Ignore the first line (the one starting with ">") and all whitespace
 *  Treat every other character as a nucleotide, case-insensitive
 *  For simplicity, just treat anything other than aAcCgG as a T.
 *  Again in the interest of simplicity, we are only going to look at the long chromosome,
 * @author Yejin Yoon
 *
 */
public class ORFParser {
	
	public static void parse(HashSet<Integer> A, HashSet<Integer> T,
								HashSet<Integer> C, HashSet<Integer> G) {
		
	}
}
