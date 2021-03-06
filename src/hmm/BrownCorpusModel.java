package hmm;

/**
 * Special HMM for 'Brown Corpus' text which has all parameters of lambda predetermined
 * @author Shayan Asgari, Sonnan Naeem 
 *
 */
public class BrownCorpusModel extends Model 
{
	
	public double[][] A;
	public double[][] B;
	public double[] pi;
	
	
	/**
	 * Number of observations to fetch
	 */
	public static final int LENGTH_OF_TEXT = 50000;
	
	/**
	 * Initialize the initial A, B, PI matrices to predetermined values
	 */
	public static final double[][] A_ARRAY = 
		{
			
											{0.47468, 0.52532},
											{0.51656, 0.48344}
											
		};
	public static final double[] PI = {0.51316,0.48684};
	
	public static final double[][] B_ARRAY = 
		{
											  {   0.03735, 
												  0.03408, 
												  0.03455,
												  0.03828,
												  0.03782,
												  0.03922,
												  0.03688,
												  0.03408, 
												  0.03875,
												  0.04062, 
												  0.03735, 
												  0.03968, 
												  0.03548,
												  0.03735, 
												  0.04062,
												  0.03595,
												  0.03641,
												  0.03408,
												  0.04062, 
												  0.03548, 
												  0.03922,
												  0.04062,
												  0.03455 ,
												  0.03595,
												  0.03408, 
												  0.03408,
												  0.03688
											  },
											  
											  { 
												  0.03909,
												  0.03537,
												  0.03537,
												  0.03909,
												  0.03583,
												  0.03630,
												  0.04048,
												  0.03537,
												  0.03816,
												  0.03909,
												  0.03490,
												  0.03723,
												  0.03537,
												  0.03909,
												  0.03397,
												  0.03397,
												  0.03816,
												  0.03676,
												  0.04048,
												  0.03443,
												  0.03537,
												  0.03955,
												  0.03816,
												  0.03723,
												  0.03769,
												  0.03955,
												  0.03397,												  
											  }
	};
										 
	public BrownCorpusModel()
	{
		super(A_ARRAY,B_ARRAY,PI);
	}
}
