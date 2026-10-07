#pragma once
#include <QStringList>

inline QStringList experimentFieldOrder() {
    return {"NAME", "PARAMETERS_TO_SOLVE_FOR", "REACTIONS", "SPECIES",
        "SPECIE_INLET_CONC_UNITS", "SPECIE_MODEL_CONC_UNITS", "EDGES", "WIDTH",
        "DEFAULT_NORMALIZATION", "LOW_REF_LEFT", "LOW_REF_RIGHT", "HIGH_REF_LEFT",
        "HIGH_REF_RIGHT", "CHANNEL_WIDTH", "CHANNEL_HEIGHT", "CHANNEL_LENGTH", "ENTRANCE_FLOWRATE"};
}

inline QString experimentFieldLabel(QString name, bool includeDimensionUnit = false) {
    if (name == "DEFAULT_NORMALIZATION") return "Normalization";
    const bool dimension = name.startsWith("CHANNEL_");
    name = name.toLower().replace('_', ' ');
    if (!name.isEmpty()) name[0] = name[0].toUpper();
    return name + (dimension && includeDimensionUnit ? " (m)" : "");
}
