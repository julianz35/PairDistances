/**
 * Program to allow the calculation of evolutionary distances for protein sequences.
 *
 * Created by julianzaugg on 17/06/2015.
 */
import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JTT;
import bn.ctmc.matrix.LG;
import bn.ctmc.matrix.Dayhoff;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.FastaReader;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;


public class PairDistances {

    private static SubstModel MODEL = new JTT();
    private String INPUT_FILENAME;
    private String OUTPUT_FILENAME;
    private boolean SKIP_OVERHANG = false;
    private double START_TIME = 0.0;
    private double END_TIME = 1.0;
    private int GA_ITER = 1000;
    private double GA_ALPHA = 0.0001; //Constant in gradient ascent, scales
    private double GA_DELTA = 0.01; //Step size in gradient ascent
    private double GA_GRADIENT_THRESHOLD = 0.1;
    private int SCAN_STEPS = 5;
    private boolean EXHAUSTIVE = false;

    private HashMap<Double, double[][] > rateMatrices = new HashMap();
    private int ROUNDING_MULTIPLIER = 1000;

    public PairDistances(){

    }

    public static void main(String args[]){
        PairDistances pd = new PairDistances();

//        String inputFastaFilename = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary Pathway/Data/Other/Honours_data/pos_neg_seqs_reduced.txt";
        String inputFastaFilename = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary Pathway/Data/Other/Li et al/CYP102A_chim_seqs_reduced.txt";
//        String inputFastaFilename = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary Pathway/" +
//                "Results/Experimental clustering/CASTing/ANEH_mutseqs_fasta.txt";
//        String inputFastaFilename = "/Users/julianzaugg/Documents/University/Phd/Projects/NLS/Data/Sequences/Bipartite.txt";
//        String inputFastaFilename = "/Users/julianzaugg/Documents/University/Phd/Projects/NLS/Data/Sequences/Monopartite.txt";
        String outputFilename = "/Users/julianzaugg/Downloads/CYP102A_chim_seqs_reduced_results27.txt";
        Double startTime = 0.1;
        Double endTime = 4.0;
        Integer GAIterLimit = 1000;
        Double GARate = 0.0001;
        Double stepDelta = 0.001;
        Boolean skipOverhangs = true;
        boolean exhaustive = false;
        double GAThresh = 1;
        int scanSteps = 5;
        int roundingMultiplier = 1000;

        pd.INPUT_FILENAME = inputFastaFilename;
        pd.OUTPUT_FILENAME = outputFilename;
        pd.START_TIME = startTime;
        pd.END_TIME = endTime;
        pd.GA_ITER = GAIterLimit;
        pd.GA_ALPHA = GARate;
        pd.GA_DELTA = stepDelta;
        pd.SKIP_OVERHANG = skipOverhangs;
        pd.GA_GRADIENT_THRESHOLD = GAThresh;
        pd.SCAN_STEPS = scanSteps;
        pd.EXHAUSTIVE = exhaustive;
        pd.ROUNDING_MULTIPLIER = roundingMultiplier;
        pd.MODEL = new JTT();
        pd.Run();
    }


    public static void main2(String args[]){
        //Create command line options
        Options options = new Options();
        options.addOption("i", "Input", true, "Filename for the input fasta file");
        options.addOption("o", "Output", true, "Filename for the results file");
        options.addOption("st", "StartTime", false, "Time to start search from, default 0");
        options.addOption("et", "EndTime", false, "Time to end search at, default 1");
        options.addOption("iter", "GAIter", false, "Max number of iterations for gradient ascent, default 1000");
        options.addOption("gaa", "GAAlpha", false, "Constant in gradient ascent function, default 0.01");
        options.addOption("gad", "GADelta", false, "Sstep size in gradient ascent, 0.01");
        options.addOption("so" ,"SkipOverhangs", false, "Whether or not to skip alignment configurations that create overhangs");
        options.addOption("gat", "GATheshold", false, "Threshold for ending a gradient ascent search, default 0.1");
        options.addOption("ps", "Prescan", false, "Number of steps/intervals in the prescan, default 5");
        options.addOption("exhaustive", false, "Perform an exhaustive search across all times, default false");
        options.addOption("rm", "RoundingMulti", false, "Value for rounding time values to N decimal places, where" +
                                    "N is the number of 0's in multiplier (must be [10,100,1000...etc], default 1000");
        options.addOption("model", false, "Evolutionary model to use, must be [JTT,LG,Dayhoff], default is JTT");

        PairDistances pd = new PairDistances();

        //Construct parser for user supplied arguments
        CommandLineParser parser = new DefaultParser();
        try {
            // parse the command line arguments
            CommandLine line = parser.parse( options, args);
            pd.INPUT_FILENAME = line.getOptionValue("i");
            pd.OUTPUT_FILENAME = line.getOptionValue("o");
            // validate that block-size has been set
            if (line.hasOption( "st")) {
                pd.START_TIME = Double.parseDouble(line.getOptionValue("st"));
                if (pd.START_TIME < 0) throw new ParseException("Start time must >= 0");
            }
            if (line.hasOption("et")){
                pd.END_TIME = Double.parseDouble(line.getOptionValue("et"));
                if (pd.END_TIME <= pd.START_TIME) throw new ParseException("End time must be > start time");
            }
            if (line.hasOption("iter")){
                pd.GA_ITER = Integer.parseInt(line.getOptionValue("iter"));
                if (pd.GA_ITER <= 0) throw new ParseException("Iteration number must be at >= 1");
            }
            if (line.hasOption("gaa")){
                pd.GA_ALPHA = Double.parseDouble(line.getOptionValue("gaa"));
            }
            if (line.hasOption("gadi")){
                pd.GA_DELTA = Double.parseDouble(line.getOptionValue("gadi"));
            }
            if (line.hasOption("so")){
                pd.SKIP_OVERHANG = Boolean.parseBoolean(line.getOptionValue("so"));
            }
            if (line.hasOption("gat")){
                pd.GA_GRADIENT_THRESHOLD = Double.parseDouble(line.getOptionValue("gat"));
                if (pd.GA_GRADIENT_THRESHOLD < 0) throw new ParseException("Gradient threshold must be >= 0.0");
            }
            if (line.hasOption("ps")){
                pd.SCAN_STEPS = Integer.parseInt(line.getOptionValue("ps"));
            }
            if (line.hasOption("exhaustive")){
                pd.EXHAUSTIVE = Boolean.parseBoolean(line.getOptionValue("exhaustive"));
            }
            if (line.hasOption("rm")){
                pd.ROUNDING_MULTIPLIER = Integer.parseInt(line.getOptionValue("exhaustive"));
            }
            if (line.hasOption("model")){
                String choice = line.getOptionValue("model");
                if (choice == "JTT") pd.MODEL = new JTT();
                else if (choice == "LG") pd.MODEL = new LG();
                else if (choice == "Dayhoff") pd.MODEL = new Dayhoff();
                else throw new ParseException("Unrecognised evolutionary model " +  choice);
            }
        }
        catch( ParseException exp ) {
            System.out.println( "Unexpected exception:" + exp.getMessage());
        }
    }


    /**
     * Run the program using the user supplied and/or default values.
     */
    public void Run(){
        //Pre-compute rate matrices
        this.computeMatrices(START_TIME, END_TIME, (int)END_TIME * ROUNDING_MULTIPLIER);
        //Load Sequences
        EnumSeq[] seqs = getFastaSeqs(INPUT_FILENAME);
        try {
            File file = new File(OUTPUT_FILENAME);
            FileWriter writer = new FileWriter(file);
            //Write header to output
            writer.write("Seq1\tSeq2\tTime\tLogLikelihood\tConfiguration\tS1Length\tS2Length\tS1Sequence\tS2Sequence\n");
            for (int i = 27; i < seqs.length; i++) {
                System.out.println(i);
                EnumSeq s1 = seqs[i];
                for (int j = 0; j <= i; j++) {
                    Object[] result = new Object[]{0.0,0, -1000000};
                    EnumSeq s2 = seqs[j];
                    //Determine which sequence is longer/shorter (just for consistency)
                    EnumSeq longerSeq = s1.length() != s2.length() ? (s1.length() > s2.length() ? s1 : s2) : s1;
                    EnumSeq shorterSeq = s1.length() != s2.length() ? (s1.length() < s2.length() ? s1 : s2) : s2;
                    //If sequences are identical, we assume optimal time is 0.0 and configuration is length - 1
                    if (s1.toString().equals(s2.toString())) {
                        int config = longerSeq.length() - 1;
                        double score = getSeqPairLL(longerSeq, shorterSeq, 0.0, config);
                        result = new Object[]{0.0, config, score};
                    }
                    else {
                        //************************
//                        for (double t = 0.0; t <= 4; t += (4-0.0)/10){
//                            for (int c =0; c < possibleConfigs(s1,s2); c++){
//                                double temp = this.getSeqPairLL(longerSeq, shorterSeq, c, t);
//                                String out = t + "\t" + c + "\t" + temp + "\t" + s1.getName() +"\t" + s2.getName() +
//                                        "\t" + s1.length() + "\t" + s2.length();
//                                System.out.println(out);
////                                writer.write(out);
//                            }
//                        }
                        //robust scan
                        Set<Integer> likelyConfigs = identifyOptimalConfigs(longerSeq, shorterSeq);
                        Object[] bestResult = new Object[]{0.0, 0, -1000000000.0}; // (time, config, score)
                        if (EXHAUSTIVE){ //if performing an exhaustive search across all pre-calculated time points
                            for (double t : rateMatrices.keySet()) {
                                for (int oc : likelyConfigs) {
                                    double score = getSeqPairLL(longerSeq, shorterSeq, t, oc);
                                    if (score > (double)bestResult[2]) {
                                        bestResult[0] = t;
                                        bestResult[1] = oc;
                                        bestResult[2] = score;
                                    }
                                }
                            }
                        }
                        else{ //if performing a gradient ascent optimised search
                            //Find the best configuration at each pre-scan time interval
                            for (int oc : likelyConfigs){
                                Object[] curResult = gradientAscent(longerSeq,shorterSeq, START_TIME, oc, GA_ITER);
//                                double time = (double)curResult[0];
//                                int config = (int)curResult[1];
                                double score = (double)curResult[2];
                                if (score > (double)bestResult[2]){
                                    bestResult = curResult;
                                }
                            }
                        }
                        result = bestResult;
                    }
                    //Construct output string
                    String result_string = s1.getName() + "\t" + s2.getName() + "\t" + result[0] + "\t" + result[2] +
                            "\t" + result[1] + "\t" + s1.length() + "\t" + s2.length() +
                            "\t" + s1.toString() + "\t" + s2.toString() + "\n";
                    writer.write(result_string); //Write output to file
                    writer.flush(); //flush the writer to keep file up to date (to avoid partial writes)
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Identify the alignment configurations for a pair of sequences, at specified time intervals,
     * that have the best score. Number of and value of intervals are dictated by {@link PairDistances#START_TIME},
     * {@link PairDistances#END_TIME} and {@link PairDistances#SCAN_STEPS}.
     * @param longerSeq
     * @param shorterSeq
     * @return Set of highest scoring configuration indexes (Integers)
     */
    public Set<Integer> identifyOptimalConfigs(EnumSeq longerSeq, EnumSeq shorterSeq){
        Set bestConfigs = new HashSet();
        int nConfigs = possibleConfigs(longerSeq, shorterSeq);
        double stepSize = (END_TIME - START_TIME)/SCAN_STEPS;
        int[] drange = getInDomainConfigRange(longerSeq,shorterSeq);
        for (double time = START_TIME; time <= END_TIME; time += stepSize) {
            int bestConfig = 0;
            double bestConfigScore = -100000000;
            for (int config = 0; config < nConfigs; config ++) {
                //If skipping overhang configurations, only consider relevant in-domain configurations
                if (SKIP_OVERHANG)
                    if (config < drange[0] || config > drange[1])
                        continue;
                double score = this.getSeqPairLL(longerSeq, shorterSeq, time, config);
                if (score > bestConfigScore){
                    bestConfig = config;
                    bestConfigScore = score;
                }
            }
            bestConfigs.add(bestConfig);
        }
        return bestConfigs;
    }

    /**
     * Returns the total number of alignment configurations possible for a pair of sequences.
     * @param longerSeq
     * @param shorterSeq
     * @return The number of alignments possible
     */
    public int possibleConfigs(EnumSeq longerSeq, EnumSeq shorterSeq){
        int overhang = shorterSeq.length() - 1;
        int variation = longerSeq.length() - shorterSeq.length();
        return variation + overhang * 2;
    }

    /**
     * Returns the minimum and maximum configuration indexes for where the shorter sequence will be aligned
     * within the domain of the longer sequence, i.e., no overhang.
     * @param longerSeq
     * @param shorterSeq
     * @return array of {int, int}, where array[0] <= array[1]
     */
    public int[] getInDomainConfigRange(EnumSeq longerSeq, EnumSeq shorterSeq){
        int min = shorterSeq.length() - 1;
        int max  = longerSeq.length() - 1;
        return new int[]{min,max};
    }

    /**
     * Performs a gradient ascent search to find the time point at which two sequences are most likely
     * to share a common ancestor at the user specified alignment configuration. Is considered converged when the
     * gradient value is less than the user specified threshold {@link PairDistances#GA_GRADIENT_THRESHOLD}.
     * @param longerSeq
     * @param shorterSeq
     * @param time A time value to search from
     * @param config An alignment configuration index
     * @param iterations Number of iterations left for performing search
     * @return An object array, where [0] = the optimal time value, [1] the supplied configuration value and
     * [2] the optimal log likelihood cost
     */
    public Object[] gradientAscent(EnumSeq longerSeq, EnumSeq shorterSeq, double time, int config, int iterations){
        double cost = getSeqPairLL(longerSeq, shorterSeq, time, config); //cost at current time
        double nextCost_t = getSeqPairLL(longerSeq, shorterSeq, time + GA_DELTA, config); //cost at time + delta
        double gradient_t = (nextCost_t - cost)/GA_DELTA;
        double temp_theta_t = time + GA_ALPHA * gradient_t; //calculate next time value
        double delta_temp = GA_DELTA;
        //If the next time goes below 0, we continually half the step size and recalculate the cost and next time
        while (temp_theta_t < 0){
            delta_temp /= 2;
            nextCost_t = getSeqPairLL(longerSeq, shorterSeq, time + delta_temp, config); //cost at time + delta
            gradient_t = (nextCost_t - cost)/delta_temp;
            temp_theta_t = time + GA_ALPHA * gradient_t;
        }
//        System.out.println(longerSeq.getName() + "\t" +
//                            shorterSeq.getName() + "\t" +
//                            cost + "\t" +
//                            nextCost_t + "\t" +
//                            gradient_t + "\t" +
//                            time + "\t" +
//                            temp_theta_t + "\t" +
//                            iterations);
        if(Math.abs(gradient_t) < GA_GRADIENT_THRESHOLD  || iterations == 0){
//            double time_orig = reverseMeanNormalise(time, START_TIME, END_TIME);
            return new Object[]{time, config, cost};
        }
        else {
            double nextTime = temp_theta_t;
            return gradientAscent(longerSeq, shorterSeq, nextTime, config, iterations - 1);
        }
    }

    @Deprecated
    public double GASeqPairLL(EnumSeq longerSeq, EnumSeq shorterSeq, double time_normed, int config){
        double time_orig = reverseMeanNormalise(time_normed, START_TIME, END_TIME);
        //For each configuration (assumed already
        return getSeqPairLL(longerSeq, shorterSeq, time_orig, config);
    }

    /**
     * Returns the log-likelihood score for the alignment of a pair of sequences at a specified alignment configuration
     * and time. Used primarily in {@see PairDistances#gradientAscent} as a cost function.
     * @param longerSeq
     * @param shorterSeq
     * @param configuration
     * @param time
     * @return A log-likelihood score (double)
     */
    public double getSeqPairLL(EnumSeq longerSeq, EnumSeq shorterSeq, double time, int configuration){

        //Ensure time is rounded appropriately and rate matrix has been calculated
        double roundedTime = (double)Math.round(time * ROUNDING_MULTIPLIER)/(double)ROUNDING_MULTIPLIER;
        if (!rateMatrices.containsKey(roundedTime)){
            computeMatrices(0, roundedTime, 1);
        }

        int longLength = longerSeq.length();
        int shortLength = shorterSeq.length();
        int overhang = shortLength - 1;
        double sumTotal = 0.0;

        //shorter sequence overhangs on left
        if (configuration < overhang) {
            for (int position = 0; position < longLength + overhang - configuration; position++) {
                //gap in longer
                if (position < overhang - configuration) {
                    sumTotal += this.getLogLikelihood(null, shorterSeq.get()[position], roundedTime);
                }
                //content in both
                else if (position >= overhang - configuration && position <= overhang) {
                    sumTotal += this.getLogLikelihood(longerSeq.get()[position - (overhang - configuration)],
                            shorterSeq.get()[position], roundedTime);
                }
                //gap in shorter
                else if (position > overhang) {
                    sumTotal += this.getLogLikelihood(longerSeq.get()[position - (overhang - configuration)], null,
                                                                                                       roundedTime);
                }
            }
        }
        //shorter within domain of longer sequence
        else if (configuration >= overhang && configuration < longLength) {
            for (int position = 0; position < longLength; position++) {
                //gap in shorter
                if (position < configuration - overhang || position > configuration) {
                    sumTotal += this.getLogLikelihood(longerSeq.get()[position], null, roundedTime);
                }
                //content in both
                else if (position >= configuration - overhang && position <= configuration) {
                    sumTotal += this.getLogLikelihood(longerSeq.get()[position],
                            shorterSeq.get()[position - (configuration - overhang)], roundedTime);
                }
            }
        }
        //shorter sequence overhangs on the right
        else if (configuration >= longLength - overhang) {
            for (int position = 0; position <= configuration; position++) {
                //gap in shorter
                if (position < longLength - (longLength + overhang) + configuration) {
                    sumTotal += this.getLogLikelihood(longerSeq.get()[position], null, roundedTime);
                }
                //gap in longer
                else if (position >= (longLength - (longLength + overhang) + configuration) && position < longLength){
                    sumTotal += this.getLogLikelihood(longerSeq.get()[position], shorterSeq.get()[position -
                                                                              (configuration - overhang)], roundedTime);
                }
                //content in both
                else if (position >= longLength) {
                    sumTotal += this.getLogLikelihood(null, shorterSeq.get()[position - (configuration - overhang)],
                                                                                                        roundedTime);
                }
            }
        }
        return sumTotal;
    }

    /**
     * Compute and store conditional probability matrices using evolutionary model - {@link PairDistances#MODEL}
     * @param startTime Start of time range
     * @param endTime End of time range
     * @param steps Steps in range
     */
    public void computeMatrices(double startTime, double endTime, int steps){
        double stepSize = (endTime-startTime)/steps;
        if (stepSize == 0) stepSize = 1;
        for (double i = startTime; i <= endTime; i += stepSize){
            double roundedTime = (double)Math.round(i * ROUNDING_MULTIPLIER) / (double)ROUNDING_MULTIPLIER;
            if (!rateMatrices.containsKey(roundedTime))
                rateMatrices.put(roundedTime, MODEL.getProbs(roundedTime));
        }
    }

    @Deprecated
    /**
     * Return mean normalised version of a value (in our case, the mean is just the (max-min/2)
     * @param value
     * @param rangeMin
     * @param rangeMax
     * @return
     */
    public double meanNormalise(double value, double rangeMin, double rangeMax){
        double range = rangeMax - rangeMin;
//        return value/(rangeMax);
        return (value - range/2)/(range);
    }

    @Deprecated
    public double reverseMeanNormalise(double value, double rangeMin, double rangeMax){
        double range = rangeMax - rangeMin;
//        return value * rangeMax;
        return value * range + range/2;
    }

    /**
     * Get the log-likelihood score for the probability of two aligned amino acids sharing the same ancestor.
     * @param aminoacid1
     * @param aminoacid2
     * @param time
     * @return Log-likelihood score representing probability
     */
    public double getLogLikelihood(Object aminoacid1, Object aminoacid2, Double time) {
        double result = 0.0;
        Object[] values = MODEL.getDomain().getValues();
        if (aminoacid1 == null && aminoacid2 == null) {
            for (Object parent : values) {
                double p_par = MODEL.getProb(parent);
                for (Object aa1 : values) {
                    double p_aa1 = MODEL.getProb(aa1, parent, rateMatrices.get(time));
                    for (Object aa2 : values) {
                        double p_aa2 = MODEL.getProb(aa2, parent, rateMatrices.get(time));
                        if (aa1 ==parent || aa2 == parent) continue;
                        result += p_par * p_aa1 * p_aa2;
                    }
                }
            }
        }
        else if (aminoacid1 != null && aminoacid2 != null) {
            for (Object parent : values) {
                double p_par = MODEL.getProb(parent);
                double p_aa1 = MODEL.getProb(aminoacid1, parent, rateMatrices.get(time));
                double p_aa2 = MODEL.getProb(aminoacid2, parent, rateMatrices.get(time));
                result += p_par * p_aa1 * p_aa2;
            }
        } else if (aminoacid1 == null && aminoacid2 != null) {
            for (Object parent : values) {
                double p_par = MODEL.getProb(parent);
                double p_aa2 = MODEL.getProb(aminoacid2, parent, rateMatrices.get(time));
                for (Object aa1 : values) {
                    if (aa1 == parent) continue;
                    double p_aa1 = MODEL.getProb(aa1, parent, rateMatrices.get(time));
                    result += p_par * p_aa1 * p_aa2;
                }
            }
        } else if (aminoacid1 != null && aminoacid2 == null) {
            for (Object parent : values) {
                double p_par = MODEL.getProb(parent);
                double p_aa1 = MODEL.getProb(aminoacid1, parent, rateMatrices.get(time));
                for (Object aa2 : values) {
                    if (aa2 == parent) continue;
                    double p_aa2 = MODEL.getProb(aa2, parent, rateMatrices.get(time));
                    result += p_par * p_aa1 * p_aa2;
                }
            }
        }
        return Math.log(result);
    }


    /**
     * Load protein sequences from a fasta file.
     * @param filename
     * @return List of enumerable protein sequences
     */
    public static EnumSeq[] getFastaSeqs(String filename){
        try {
            FastaReader freader = new FastaReader(filename, Enumerable.aacid);
            EnumSeq[] seqs = freader.load();
            return seqs;
        } catch (IOException e){
            e.printStackTrace();
            return null;
        }
    }

}
