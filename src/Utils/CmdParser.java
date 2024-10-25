// src/Utils/CmdParser.java
package Utils;

import java.util.HashMap;
import java.util.Map;

public class CmdParser {
    private final Map<String, String> arguments = new HashMap<>();

    public void declareArgument(String argName) {
        arguments.put(argName, null);
    }

    public void parse(String[] args) {
        for (int i = 0; i < args.length; i++) {
            if (arguments.containsKey(args[i]) && i + 1 < args.length) {
                arguments.put(args[i], args[i + 1]);
                i++; // Skip the value
            }
        }
    }

    public String getArgumentValue(String argName) {
        return arguments.get(argName);
    }
}