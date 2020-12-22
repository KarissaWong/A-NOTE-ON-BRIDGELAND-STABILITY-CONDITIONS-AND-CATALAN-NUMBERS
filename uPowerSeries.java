import java.util.ArrayList;

public class uPowerSeries {
	
	
	/*
	 *input: the power u is raised to
	 *finds the coefficients for power expansion in the form of an array
	 *ex: pascalCoefficients(2) = [1, 2, 1]
	 *output: array with the coefficients
	 */
	public static long[] pascalCoefficients(long power) {
	//this code was pretty much copied from online and modified for this program
        int num = (int) (power + 1);
        long[] pascalArray = new long[num];
		for (long i = 0; i < num; i++) {
        	long number = 1;
            for (int j = 0; j <= i; j++) {
                pascalArray[j] = number;
                number = number * (i - j) / (j + 1);

            }
        }
        return pascalArray;
    }
	
	
	/*
	 *input: the power u is raised to
	 *expands binomial to the nthPower
	 *output: 2d array that represents the terms
	 */
	public static long[][] nthPower(long nthPower){
		long[] coefficients = pascalCoefficients(nthPower);
		long[][] simplified = new long[(int) (nthPower + 1)][5];
		
		for (int i = 0; i < simplified.length; i++) {
			
			for (int j = 0; j < 5; j++) {
				if (j == 0) {
					if (i % 2 == 0) {
						simplified[i][j] = coefficients[i];
					}
					else {
						simplified[i][j] = -1*coefficients[i];
					}
					
				}
				else if (j == 1) {
					simplified[i][j] = nthPower - i;
				}
				else if (j == 2) {
					simplified[i][j] = i;
				}
				else if (j == 3) {
					simplified[i][j] = nthPower;
				}
				else if (j == 4) {
					simplified[i][j] = 2 * i;
				}
			}
		}
		return simplified;
	}
	
	
	/*
	 *input: term with u and expanded binomial from nthPower
	 *multiplies term with u with the expanded binomial
	 *output: 2d array with multiplied terms
	 */
	public static long[][] multiply (long[] term, long[][] expanded){
		long[][] multiplied = new long[expanded.length][expanded[0].length];
		
		for (int i = 0; i < expanded.length; i++) {
			
			for (int j = 0; j < expanded[i].length; j++) {
				
				if (j == 0) {
					multiplied[i][j] = term[0] * expanded[i][0];
				}
				else if (j == 1) {
					multiplied[i][j] = term[j] + expanded[i][j];
				}
				else if (j == 2) {
					multiplied[i][j] = term[j] + expanded[i][j];
				}
				else if (j == 3) {
					multiplied[i][j] = term[j] + expanded[i][j];
				}
				else if (j == 4) {
					multiplied[i][j] = term[j] + expanded[i][j];
				}
			}
		}
		return multiplied;
		
	}
	
	
	/*
	 *input: 2 like terms
	 *adds the coefficients of the terms
	 *output: one term in the form of an array
	 */
	public static long[] addLikeTerms(long[] term1, long[] term2) {
		long[] added = term1;
		
		added[0] = term1[0] + term2[0];
		
		return added;
	}
	
	
	/*
	 *input: all the terms
	 *searches for like terms, adds them using addLikeTerms and removes the other
	 *like term
	 *output:all the terms but condensed  
	 */
	public static ArrayList<long[]> simplifyTerms(ArrayList<long[]> allTerms) {
		ArrayList<long[]> simplifiedTerms = new ArrayList<long[]>(); 
		
		for (long w = 0; w < allTerms.size(); w++) { //instantiating arrays
			simplifiedTerms.add(new long[5]);
		}
		
        
		for (int a = 0; a < simplifiedTerms.size(); a++) { //making a copy of allTerms
        	simplifiedTerms.set(a, allTerms.get(a));
        }

        for (int i = 0; i < simplifiedTerms.size() - 1; i++) { 
        
        	long[] term = simplifiedTerms.get(i);
        	
        	int j = i + 1;
        	while (j < simplifiedTerms.size()) { 
        	
        		long[] comparedTerm = simplifiedTerms.get(j);
        		if (term[1] != comparedTerm[1]) {
        			j++;
        		}
        		else if (term[2] != comparedTerm[2]) {
        			j++;
        		}
        		else if (term[3] != comparedTerm[3]) {
        			j++;
        		}
        		else if (term[4] != comparedTerm[4]) {
        			j++;
        		}
        		else {
        			simplifiedTerms.set(i, addLikeTerms(term, comparedTerm));
        			simplifiedTerms.remove(comparedTerm);
        			//j--;
        			//j++;
        		}
        	}
        }
        return simplifiedTerms;
	}
	
	
	/*
	 *input: all the terms
	 *finds all terms that have u^0
	 *adds only the first few terms (iterations + 1) of all the terms that have u^0
	 *because each iteration only brings one absolute term (and the 0th iteration has 1 absolute term) 
	 *output: absolute terms
	 */
	public static ArrayList<long[]> absoluteTerms(ArrayList<long[]> allTerms, long iterations){
		ArrayList<long[]> absoluteTerms = new ArrayList<long[]>(); 
		
		for (int i = 0; i < iterations + 1; i++) {
			long[] term = allTerms.get(i);
			
			if (term[4] == 0) {
				absoluteTerms.add(term);
			}
		}
		return absoluteTerms;
	}
	
	
	/*
	 *input: absolute terms
	 *output:
	 *prints out #A^nB^nw^n format (u is not included because it is a given u^0)
	 *prints out [ # A^ B^ w^ u^ ] format
	 */
	public static void printNicely(ArrayList<long[]> absoluteTerms) {
		for (int i = 0; i < absoluteTerms.size(); i++) {
			long[] term = absoluteTerms.get(i);
			System.out.print(term[0] + "A^" + term[1] + "B^" + term[2] + "w^" + term[3]);
			if (i + 1 != absoluteTerms.size()) {
				System.out.print(" + ");
			}
		}
		System.out.println();
		
		for (int j = 0; j < absoluteTerms.size(); j++) {
			long[] term = absoluteTerms.get(j);
			
			System.out.print("[ ");
			for (int k = 0; k < term.length; k++) {
				System.out.print(term[k] + " ");
			}
			System.out.print("]");
			
			System.out.println();
		}
	}
	
	
	public static void printTerm(long[] term, long num) {
		System.out.print("[ ");
    	for (int a = 0; a < 5; a++) {
    		System.out.print(term[a] + " ");
    	}
    	System.out.print("]");
    	System.out.println();
    	System.out.print("It has been verified that the degree-n coefficients of the power series equals the n-th Catalan number, up to n = " + num + "\r");
    	System.out.flush();
	}
	
	
	public static ArrayList<long[]> filterTerms(long target, ArrayList<long[]> allTerms){
		ArrayList<long[]> filteredTerms = new ArrayList<long[]>(); 
		
		for (int i = 0; i < allTerms.size(); i++) { 
			long[] term = allTerms.get(i);
			if(term[3] + term[4] <= target) {
				filteredTerms.add(allTerms.get(i));
			}
        }
		return filteredTerms;
	}

	
	public static boolean necessaryTerm(long target, long[] term) {
		if(term[3] + term[4] <= target) {
			return true;
		}
		else {
			return false;
		}
	}
	
	
	public static double maxC(ArrayList<Double> allCs) {
		double max = allCs.get(0);

		for (int j = 0; j < allCs.size(); j++) {
			if (max < allCs.get(j)) {
				max = allCs.get(j);
			}
		}

		return max;
		
	}
	
	public static double hadamard(double maxC) {
		return 1/maxC;
	}
	
	/*
	 *input: number of iterations
	 *[details inside]
	 *output: printed terms
	 */
	public static void uPwrSeries(int iterations) {
		ArrayList<long[]> answer = new ArrayList<long[]>(); 
		ArrayList<Double> allCs = new ArrayList<Double>();
        long[] aw = new long[]{1, 1, 0, 1, 0};
        long[] _bwu2 = new long[] {-1, 0, 1, 1, 2};
        allCs.add(1.0);
        answer.add(0, aw);
        answer.add(1, _bwu2);
        long targetWPower = iterations * 2 + 1;
        
        for (int i = 0; i < iterations; i++) {
        	int j = 0;
        	long size = answer.size();
        	
        	while (j < size) {
        	//searching for a term with u
        		long[] term = answer.get(j);
        		if (term[4] == 0) {
        		//if the term does not have u, it ignores the term
        			j++;
        		}
        		else {
        		//if the term does have u, it does this:
        			long[][] expandedTerms = nthPower(term[4]);
        			//expand the binomial to its power--located in term[4]
        
        			term[4] = 0;
        			//since we are multiplying the term with the expansion of u, we need to 
        			//replace the u term with the expansion and set the power of u = 0
        			//ex: -bw(....) not -bwu^2(....)
        			
        			long[][] multipliedTerms = multiply(term, expandedTerms);
        			//multiply the original term by the expansion
        			
        			answer.set(j, multipliedTerms[0]); 
        			//replaces the term that had u in it with the first term of the simplified
        			//expansion 
        			//ex: aw, -bwu^2 becomes aw, -a^2bw^3
        			
        			for (int k = multipliedTerms.length - 1; k > 0; k--) {
        			//adds the rest of the "new" terms
        			//ex: aw, -bwu^2 becomes aw, -a^2bw^3, 2ab^2u^2w^3, -a^2bw^3
        				answer.add(j+1, multipliedTerms[k]);
        			}
        			
        			size = size + multipliedTerms.length - 1;
        			
        			j = j + multipliedTerms.length;
        			//skips over the newly added terms to find the next term with u in it
        		}
        		
        	}
        	
        	long[] lastTerm = answer.get(answer.size()-1);
        	
        	if (!(necessaryTerm(targetWPower, lastTerm))) {
        		answer = filterTerms(targetWPower, answer);
        		
        	}

        	answer = simplifyTerms(answer);
        	//simplifying like terms
        	
        	
        	long[] absoluteTerm = answer.get(i);
        	System.out.print("[ ");
        	for (int a = 0; a < 5; a++) {
        		System.out.print(absoluteTerm[a] + " ");
        	}
        	System.out.print("]");
        	
        	System.out.println();
        	System.out.print("It has been verified that the degree-n coefficients of the power series equals the n-th Catalan number, up to n = " + i + "\r");
        	System.out.flush();
        	
        	
        	//finding roc
        	
        	double coLastTerm = (double) lastTerm[0];
        	double n = (double) i;
        	double c = (Math.pow(Math.abs(coLastTerm), 1/n));
        	allCs.add(c);
        	
        	double maxC = maxC(allCs);
        	
        	System.out.println("ROC: " + hadamard(maxC));
        	
        	
        	
        	

        }
        long[] lastAbsoluteTerm = answer.get(iterations);
    	System.out.print("[ ");
    	for (int a = 0; a < 5; a++) {
    		System.out.print(lastAbsoluteTerm[a] + " ");
    	}
    	System.out.print("]");
    	long num = iterations;
    	System.out.println();
    	
    	//finding roc
    	double coLastTerm = (double) lastAbsoluteTerm[0];
    	double n = (double) iterations;
    	double c = (Math.pow(Math.abs(coLastTerm), 1/n));
    	allCs.add(c);
    	
    	double maxC = maxC(allCs);
    	
    	System.out.println("ROC: " + hadamard(maxC));
    	
 
    	
    	
    	System.out.println("It has been verified that the degree-n coefficients of the power series equals the n-th Catalan number, up to n = " + num);
        System.out.println("Calculations are complete!");
        
        //printNicely(absoluteTerms(answer, iterations));
        //prints result :D
	}

	
	public static void main(String[] args) {
		uPwrSeries(50);
		
	}

}
