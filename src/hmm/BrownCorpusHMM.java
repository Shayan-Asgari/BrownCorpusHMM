package hmm;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;

import javax.swing.JOptionPane;

import java.lang.Math;
import java.text.DecimalFormat;

/**
 * Class which implements and runs a HMM model from the first 50k characters of 'Brown Corpus'
 * @author Shayan Asgari, Sonnan Naeem 
 *
 */
public class BrownCorpusHMM 
{
	private int maxIterations;
	private int observationSequenceLength;
	private BrownCorpusModel model;
	
	//Used for scaling alpha values
	public double[] c;
	
	
	public BrownCorpusHMM(int maxIterations, int observationSequenceLength)
	{
		this.maxIterations = maxIterations;
		this.observationSequenceLength = observationSequenceLength;
		this.model = new BrownCorpusModel();
	}
	
	/**
	 * Helper method which transposes a 2D matrix
	 * @param matrix to transpose
	 * @return transposedMatrix 
	 */
	public static double[][] transposeMatrix(double[][] matrix)
	{
	    int m = matrix.length;
	    int n = matrix[0].length;

	    double[][] transposedMatrix = new double[n][m];

	    for(int x = 0; x < n; x++) {
	        for(int y = 0; y < m; y++) {
	            transposedMatrix[x][y] = matrix[y][x];
	        }
	    }
	    return transposedMatrix;
	}
	
	/**
	 * Get observations from first 50k characters of 'Brown Corpus' only including (a-z, ' ')
	 * @param url the url of the 'Brown Corpus' .txt file
	 * @return observations an ArrayList<Integer> containing integers representing 
	 * letters(0-26) and spaces(27)
	 */
	public int[] getObservations(String url) throws Exception
	{
		 ArrayList<Integer> observations = new ArrayList<Integer>();
		 URL oracle = new URL(url);
		 BufferedReader in = new BufferedReader(new InputStreamReader(oracle.openStream()));
		 String inputLine;
		 int charactersReadSoFar = 0;
		 while ((inputLine = in.readLine()) != null)
		 {
			 if(!inputLine.isBlank())
			 {
				 char[] arr = inputLine.toCharArray();
				 for(char c: arr)
				 {
					 char b = Character.toLowerCase(c);
					 if(Character.isLetter(b))
					 {
						 observations.add(b - 'a'); // a =0, z = 25
						 charactersReadSoFar++;
					 }
					 else if(Character.isSpaceChar(b))
					 {
						 observations.add(26); //space = 26
						 charactersReadSoFar++;
					 }
					 if(charactersReadSoFar == this.observationSequenceLength)
					 {
						 int[] o= new int[observationSequenceLength];
						 for(int z = 0; z<observationSequenceLength; z++)
						 {
							 o[z] = observations.get(z);
						 }
							 
						 return o;	
					 }
				 }
			 }
		 }
		 in.close();
		 return null;	
	}
	
	
	/**
	 * Alpha pass implementation outputs an alpha array that helps compute score of model
	 * @param A the state transition probabilities
	 * @param B the observation probability matrix
	 * @param pi the initial state distribution
	 * @param O the observation sequence (50k Brown Corpus text characters)
	 * @return alpha the array of alpha values produced from alpha pass
	 */
	public double[][] alphaPass(double[][] A, double[][] B, double[] pi, int[] O)
	{
		int T = O.length;
		int N = pi.length;
		
		double[][] alpha = new double[T][N];
		
		double[] c = new double[T];
		
		//Compute alpha[0][i]
		c[0] = 0;
		for (int i = 0; i < N; i++) {
			alpha[0][i] = pi[i] * B[i][O[0]];
			c[0] = c[0] + alpha[0][i];
		}
		
		//Scale the alpha[0][i]
		c[0] = 1.0 / c[0];
		for (int i = 0; i < N; i++) {
			alpha[0][i] = c[0] * alpha[0][i];
		}
		
		//Compute alpha[t][i]
		for (int t = 1; t < T; t++) {
			c[t] = 0;
			for (int i = 0; i < N;i++) {
				alpha[t][i] = 0;
				for (int j = 0; j < N; j++)
				{
					alpha[t][i] += alpha[t-1][j] * A[j][i];
					
				}
				// [1,27
				alpha[t][i] =  alpha[t][i] * B[i][O[t]];
				c[t] += alpha[t][i];
			}
			
			//Scale alphas[t][i]
			c[t] = 1 / c[t];
			for (int i = 0; i < N; i++) {
				alpha[t][i] = c[t] * alpha[t][i];
			}
		}
		
		this.c = c;
		return alpha;
	}
	
	/**
	 * Beta pass implementation outputs a beta array that helps compute the best number of hidden states
	 * @param A the state transition probabilities
	 * @param B the observation probability matrix
	 * @param pi the initial state distribution
	 * @param O the observation sequence (50k Brown Corpus text characters)
	 * @return alpha the array of beta values produced from beta pass
	 */
	public double[][] betaPass(double[][]A, double[][]B, double[] pi, int[] O) //TODO Do we get c[] from alpha or from beta pass
	{	
		int T = O.length;
		int N = pi.length;
		
		double[][] beta = new double[T][N];
		double[] c = this.c;
		
		//Let beta[T-1][i] = 1, scaled by c[t]
		for (int i = 0; i < N; i++) {
			beta[T-1][i] = c[T-1];
		}
		
		//Beta pass
		for (int t = T - 2; t >= 0; t--) {
			for (int i = 0; i < N; i++) {	
				beta[t][i] = 0;
				for (int j = 0; j < N; j++) {
					beta[t][i] = beta[t][i] + A[i][j] * B[j][O[t+1]] * beta[t+1][j];
				}
				//Scale beta[t][i] with same scale factor as alpha[t][i]
				beta[t][i] = c[t] * beta[t][i];
			}
		}
		return beta;
	}
	
	/**
	 * Build a Hidden Markov Model by using alphaPass and betaPass to help re-estimate the model
	 * @return model the final HMM model with a modified A, B, and PI matrices
	 */
	public BrownCorpusModel buildHMM() throws Exception {
		
		double[][] A = model.getA();
		double[][] B = model.getB();
		double[] pi = model.getPi();
		int[] O = this.getObservations("http://www.sls.hawaii.edu/bley-vroman/brown_nolines.txt");
		
		int N = A.length;
		int M = B[0].length;
		int T = O.length;
		
		double[][] alpha = new double[T][N];
		double[][] beta = new double[T][N];
		
		int iterations = 0;
		double oldLogProb = Double.NEGATIVE_INFINITY;
		boolean iterate = true;
		
		while (iterate) {
		
			//TODO will this work? Cause beta calls alpha, so should I call alpha? FIXED
			alpha = alphaPass(A, B, pi, O);
			beta = betaPass(A, B, pi, O);
			
			//Compute digamma and gamma
			double[][][] digamma = new double[T][N][N];
			double[][] gamma = new double[T][N];
			
			//No need to normalize digamma[t][i][j] since using scaled alpha and beta
			for (int t = 0; t < T - 1; t++) {
				for (int i = 0; i < N; i++) {
					gamma[t][i] = 0;
					for (int j = 0; j < N; j++) {
						double z = B[j][O[t+1]];
						digamma[t][i][j] = (alpha[t][i] * A[i][j] * B[j][O[t+1]] * beta[t+1][j]);
						gamma[t][i] = gamma[t][i] + digamma[t][i][j];
					}
				}
			}
			
			
			//Special case for gamma[T-1][i] (as above, no need to normalize)
			for (int i = 0; i < N; i++) {
				gamma[T-1][i] = alpha[T-1][i];
			}
			
			/* Re-estimate A, B, and pi */
			
			//Re-estimate pi
			for (int i = 0; i < N; i++) {
				pi[i] = gamma[0][i];
			}
			
			//Re-estimate A
			for (int i = 0; i < N; i++) {
				double denom = 0;
				for (int t = 0; t < T - 1; t++) {
					denom += gamma[t][i];
				}
				for (int j = 0; j < N; j++) {
					double numer = 0;
					for (int t = 0; t < T - 1; t++) {
						numer = numer + digamma[t][i][j];
					}
					A[i][j] = numer/denom;
				}
			}

			//Re-estimate B
			for (int i = 0; i < N; i++) {
				double denom = 0;
				for (int t = 0; t < T; t++) {
					denom += gamma[t][i];
				}
				for (int j = 0; j < M; j++) {
					double numer = 0;
				
					for (int t = 0; t < T; t++) {
						
						if (O[t] == j) {
							numer = numer + gamma[t][i];
						}
					}
					B[i][j] = numer/denom;
				}
			}
			
			//Compute log[P(O|gamma)]
			
			double logProb  = 0;
			for (int i = 0; i < T; i++) {
				logProb += Math.log(c[i]);
			}
			logProb *= -1;
			
			/* To iterate or not to iterate, that is the question... */
			iterations += 1;
			
			if ((iterations < this.maxIterations) && (logProb > oldLogProb)) {
				oldLogProb = logProb;
			} else {
				iterate = false;
				model.setA(A);
				model.setB(B);
				model.setPi(pi);
				System.out.println("Final Log Probability: log[P(O|Lambda)] " + logProb + "\n");
			}
		}
		return model;
	}
	
	public static void main(String[] args)
	{
		char[] lettersToRepresent = "abcdefghijklmnopqrstuvwxyz*".toCharArray();
		int indexOfLetter = 0;
		try
		{
			BrownCorpusHMM hmm = new BrownCorpusHMM(100, BrownCorpusModel.LENGTH_OF_TEXT);
			BrownCorpusModel model = hmm.buildHMM();
			System.out.println("FINAL PI MATRIX:");
			System.out.println(Arrays.toString(model.getPi()) + "\n");
			
			double[][] A = model.getA();
			System.out.println("FINAL A MATRIX: \n");
			for (int row = 0; row < A.length; row++) {
		        for (int col = 0; col < A[row].length; col++) {
		            System.out.printf("%.5f  ", A[row][col]);
		        }
		        System.out.println();
		    }
			
			System.out.println("\nFINAL B MATRIX: \n");
			double[][] B = BrownCorpusHMM.transposeMatrix(model.getB());	
			System.out.println("     1        2");
			for (int row = 0; row < B.length; row++) {
				System.out.print(lettersToRepresent[indexOfLetter++] + " ");
		        for (int col = 0; col < B[row].length; col++) {
		            System.out.printf("%.5"
		            		+ "f  ", B[row][col]);
		        }
		        System.out.println();
		    }
			
		} 
		catch (Exception e) 
		{
			JOptionPane.showMessageDialog(null, e.getMessage());
		}
	}
}
