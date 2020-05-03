package hmm;

/**
 * Generic Model class for building an HMM which initializes all of lambda's parameters
 * @author Shayan Asgari, Sonnan Naeem 
 *
 */
public class Model {
	public double[][] A;
	public double[][] B;
	public double[] pi;
	
	/**
	 * Generic Model initializer for HMM
	 * @param A
	 * @param B
	 * @param pi
	 */
	public Model (double[][] A, double[][] B, double[] pi) 
	{
		this.A = A;
		this.B = B;
		this.pi = pi;
	}

	/**
	 * Retrieve A matrix
	 * @return A
	 */
	public double[][] getA() 
	{
		return A;
	}

	/**
	 * Retrieve B matrix
	 * @return B
	 */
	public double[][] getB() {
		return B;
	}
	
	/**
	 * Retrieve PI matrix
	 * @return pi
	 */
	public double[] getPi() {
		return pi;
	}
	
	/**
	 * Set the A matrix to a new one
	 * @param a the new matrix 
	 */
	public void setA(double[][] a) 
	{
		A = a;
	}

	/**
	 * Set the B matrix to a new one
	 * @param b the new matrix 
	 */
	public void setB(double[][] b) 
	{
		B = b;
	}
	
	/**
	 * Set the PI matrix to a new one
	 * @param pi the new matrix 
	 */
	public void setPi(double[] pi) 
	{
		this.pi = pi;
	}
}
