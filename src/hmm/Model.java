package hmm;

public class Model {
	public double[][] A;
	public double[][] B;
	public double[] pi;
	
	public Model (double[][] A, double[][] B, double[] pi) {
		this.A = A;
		this.B = B;
		this.pi = pi;
	}

	/**
	 * @return the a
	 */
	public double[][] getA() {
		return A;
	}

	/**
	 * @param a the a to set
	 */
	public void setA(double[][] a) {
		A = a;
	}

	/**
	 * @return the b
	 */
	public double[][] getB() {
		return B;
	}

	/**
	 * @param b the b to set
	 */
	public void setB(double[][] b) {
		B = b;
	}

	/**
	 * @return the pi
	 */
	public double[] getPi() {
		return pi;
	}

	/**
	 * @param pi the pi to set
	 */
	public void setPi(double[] pi) {
		this.pi = pi;
	}
}
