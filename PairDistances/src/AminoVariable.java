import dat.EnumVariable;
import dat.Enumerable;
import dat.Variable;

/**
 * Extension for EnumVariable. Allows for significantly faster calculations when calling getIndex (which can
 * occur a lot).
 * Use if wishing to create a variable representing an amino acid.
 *
 * Currently unused for {@see PairDistances} as of 22/06/2015
 *
 * Created by julianzaugg on 20/01/2015.
 */
public class AminoVariable extends EnumVariable {


    public AminoVariable(Enumerable domain) {
        super(domain);
    }

    public AminoVariable(Enumerable domain, String name) {
        super(domain, name);
    }


    @Override
    public int getIndex(Object value) {
        Object[] aaValues = super.getDomain().getValues();
        for (int i =0; i < aaValues.length; i++){
            if (aaValues[i].equals(value)){
                return i;
            }
        }
        throw new RuntimeException("No value found");
    }

}
