package hmm;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;

public class ReadFromFile {
	public static final int LENGTH_OF_TEXT = 50000;
	public static ArrayList<Integer> getObservations(String url) throws Exception
	{
		 ArrayList<Integer> observations = new ArrayList<Integer>();
		 URL oracle = new URL(url);
		 BufferedReader in = new BufferedReader(new InputStreamReader(oracle.openStream()));
		 String inputLine;
		 int charactersReadSoFar = 0;
		 int inputNumber =0;
		 while ((inputLine = in.readLine()) != null)
		 {
	
			 if(!inputLine.isBlank())
			 {
				 System.out.println(inputLine);
				 char[] arr = inputLine.toCharArray();
				 for(char c: arr)
				 {
					 char b = Character.toLowerCase(c);
		
					 if(Character.isLetter(b))
					 {
						 System.out.println(b);
						 observations.add(b - 'a'); // a =0, z = 26
						 charactersReadSoFar++;
					 }
					 else if(Character.isSpaceChar(b))
					 {
						 observations.add(26); //space = 27
						 charactersReadSoFar++;
					 }
					 if(charactersReadSoFar == LENGTH_OF_TEXT)
					 {
						 System.out.println(observations.size());
						 return observations;	
					 }
				 }
			 }
			 inputNumber++;
			 if(inputNumber ==1)
			 {
				 return observations;
			 }
		 }
		 in.close();
		 return observations;		
	}
	
	public static void main(String[] args)
	{
			System.out.println('h' -'a');
		try 
		{
			ArrayList<Integer> observations = ReadFromFile.getObservations("http://www.sls.hawaii.edu/bley-vroman/brown_nolines.txt");
		System.out.println(observations);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
