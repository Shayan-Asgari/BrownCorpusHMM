package hmm;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;

import javax.swing.JOptionPane;

import java.lang.Math;
import java.text.DecimalFormat;

/* THERE ARE COMMENTS PRINTED IN THE OUTPUT OF THE PROGRAM */

/**
 * Class which implements and runs a HMM model from the first 50k characters of
 * 'Brown Corpus'
 *
 * @author Shayan Asgari, Sonnan Naeem
 *
 */
public class HMM {
	private int maxIterations;
	private HmmModel model;

	//Used for scaling alpha values
	public double[] c;

	public HMM(HmmModel model, int maxIterations) {
		this.maxIterations = maxIterations;
		this.model = model;
	}

	/**
	 * Helper method which transposes a 2D matrix
	 * 
	 * @param matrix to transpose
	 * @return transposedMatrix
	 */
	public static double[][] transposeMatrix(double[][] matrix) {
		int m = matrix.length;
		int n = matrix[0].length;

		double[][] transposedMatrix = new double[n][m];

		for (int x = 0; x < n; x++) {
			for (int y = 0; y < m; y++) {
				transposedMatrix[x][y] = matrix[y][x];
			}
		}
		return transposedMatrix;
	}

	/**
	 * Alpha pass implementation outputs an alpha array that helps compute score of
	 * model
	 * 
	 * @param A  the state transition probabilities
	 * @param B  the observation probability matrix
	 * @param pi the initial state distribution
	 * @param O  the observation sequence (50k Brown Corpus text characters)
	 * @return alpha the array of alpha values produced from alpha pass
	 */
	public double[][] alphaPass(double[][] A, double[][] B, double[] pi, int[] O) {
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
			for (int i = 0; i < N; i++) {
				alpha[t][i] = 0;
				for (int j = 0; j < N; j++) {
					alpha[t][i] += alpha[t - 1][j] * A[j][i];

				}
				alpha[t][i] = alpha[t][i] * B[i][O[t]];
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
	 * Beta pass implementation outputs a beta array that helps compute the best
	 * number of hidden states
	 * 
	 * @param A  the state transition probabilities
	 * @param B  the observation probability matrix
	 * @param pi the initial state distribution
	 * @param O  the observation sequence (50k Brown Corpus text characters)
	 * @return alpha the array of beta values produced from beta pass
	 */
	public double[][] betaPass(double[][] A, double[][] B, double[] pi, int[] O) // TODO Do we get c[] from alpha or
																					// from beta pass
	{
		int T = O.length;
		int N = pi.length;

		double[][] beta = new double[T][N];
		double[] c = this.c;

		//Let beta[T-1][i] = 1, scaled by c[t]
		for (int i = 0; i < N; i++) {
			beta[T - 1][i] = c[T - 1];
		}

		//Beta pass
		for (int t = T - 2; t >= 0; t--) {
			for (int i = 0; i < N; i++) {
				beta[t][i] = 0;
				for (int j = 0; j < N; j++) {
					beta[t][i] = beta[t][i] + A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
				}
				//Scale beta[t][i] with same scale factor as alpha[t][i]
				beta[t][i] = c[t] * beta[t][i];
			}
		}
		return beta;
	}

	/**
	 * Build/Train a Hidden Markov Model by using alphaPass and betaPass
	 * 
	 * @return model the final HMM model with a modified A, B, and PI matrices
	 */
	public HmmModel buildHMM(){

		double[][] A = model.getA();
		double[][] B = model.getB();
		double[] pi = model.getPi();

		int N = A.length;
		int M = B[0].length;
		int[] O = model.getTrainObservationSequence();
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
						double z = B[j][O[t + 1]];
						digamma[t][i][j] = (alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j]);
						gamma[t][i] = gamma[t][i] + digamma[t][i][j];
					}
				}
			}

			//Special case for gamma[T-1][i] (as above, no need to normalize)
			for (int i = 0; i < N; i++) {
				gamma[T - 1][i] = alpha[T - 1][i];
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
					A[i][j] = numer / denom;
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
					B[i][j] = numer / denom;
				}
			}

			//Compute log[P(O|gamma)]

			double logProb = 0;
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
	
		
		DataProcessor processor = new DataProcessor(
				"//Users//jay//185C//Malware-Opcodes-HMM-SVM//src//hmm//Sorted-Malicia-Opcodes//balanced//cleaman", 100,
				1);

		int M = processor.getM();
		int[] trainObservationSequence = processor.getObservationSequence();
		int[] testObservationSequence = processor.getObservationSequenceTest();
		int indexOfSymbol = 0;
		int[] symbols = new int[M];
		for(int i = 0; i<M;i++)
		{
			symbols[i] = i;
		}
		HmmModel model = new HmmModel(2, M, trainObservationSequence, testObservationSequence);
		System.out.println("******************BEFORE*****************");
		
		double[][] A = model.getA();
		System.out.println("INITIAL A MATRIX: \n");
		for (int row = 0; row < A.length; row++) {
	        for (int col = 0; col < A[row].length; col++) {
	            System.out.printf("%.5f  ", A[row][col]);
	        }
	        System.out.println();
		}
		
		System.out.println("\nINITIAL B MATRIX: \n");
		double[][] B = HMM.transposeMatrix(model.getB());	
		System.out.println("     1        2");
		for (int row = 0; row < B.length; row++) {
			System.out.print(symbols[indexOfSymbol++] + " ");
	        for (int col = 0; col < B[row].length; col++) {
	            System.out.printf("%.5"
	            		+ "f  ", B[row][col]);
	        }
	        System.out.println();
	    }
		
		//***************************************************************
		HMM cleamanModel = new HMM(model, 100);
		cleamanModel.buildHMM();
		System.out.println("******************After*****************");
	    A = model.getA();
		System.out.println("FINAL A MATRIX: \n");
		for (int row = 0; row < A.length; row++) {
	        for (int col = 0; col < A[row].length; col++) {
	            System.out.printf("%.5f  ", A[row][col]);
	        }
	        System.out.println();
		}
	    
		System.out.println("\nFINAL B MATRIX: \n");
		B = HMM.transposeMatrix(model.getB());	
		System.out.println("     1        2");
		for (int row = 0; row < B.length; row++) {
			System.out.print(symbols[indexOfSymbol++] + " ");
	        for (int col = 0; col < B[row].length; col++) {
	            System.out.printf("%.5"
	            		+ "f  ", B[row][col]);
	        }
	        System.out.println();
	    }
	}


}
