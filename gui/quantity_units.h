#pragma once
#include <QComboBox>
#include <QLineEdit>
#include <QStringList>
#include <cmath>
#include <stdexcept>

inline QString convertQuantityList(const QString& text, double factor) {
    QStringList result;
    for (const auto& token : text.split(' ', Qt::SkipEmptyParts)) {
        bool ok = false;
        const double number = token.toDouble(&ok), converted = number * factor;
        if (!ok || !std::isfinite(number) || !std::isfinite(converted)
            || (number != 0 && converted == 0))
            throw std::invalid_argument("Enter finite values within the selected units' numeric range.");
        result << QString::number(converted, 'g', 17);
    }
    return result.join(' ');
}

inline QComboBox* quantityUnits(const QString& field, QLineEdit* line) {
    auto* combo = new QComboBox;
    combo->setObjectName(field + "Units");
    if (field == "ENTRANCE_FLOWRATE") {
        combo->addItem("m³/s", 1.0); combo->addItem("L/s", 1e-3);
        combo->addItem("mL/s", 1e-6); combo->addItem("µL/s", 1e-9);
        combo->addItem("mL/min", 1e-6 / 60); combo->addItem("µL/min", 1e-9 / 60);
        combo->addItem("nL/min", 1e-12 / 60);
    } else {
        combo->addItem("m", 1.0); combo->addItem("cm", 1e-2);
        combo->addItem("mm", 1e-3); combo->addItem("µm", 1e-6); combo->addItem("nm", 1e-9);
    }
    combo->setProperty("previousIndex", 0);
    combo->setProperty("baseText", line->text());
    combo->setProperty("displayText", line->text());
    QObject::connect(combo, &QComboBox::currentIndexChanged, line, [combo, line](int index) {
        const int previous = combo->property("previousIndex").toInt();
        try {
            const auto base = line->text() == combo->property("displayText").toString()
                ? combo->property("baseText").toString()
                : convertQuantityList(line->text(), combo->itemData(previous).toDouble());
            const auto converted = convertQuantityList(base, 1 / combo->itemData(index).toDouble());
            line->setText(converted); combo->setProperty("previousIndex", index);
            combo->setProperty("baseText", base); combo->setProperty("displayText", converted);
        } catch (const std::exception& e) {
            combo->blockSignals(true); combo->setCurrentIndex(previous); combo->blockSignals(false);
            line->setToolTip(QString::fromUtf8(e.what()));
        }
    });
    return combo;
}

inline QString quantityStoredText(QLineEdit* line, QComboBox* combo) {
    if (line->text() == combo->property("displayText").toString()) {
        const auto base = combo->property("baseText").toString();
        convertQuantityList(base, 1); // Validate while retaining the original precision and formatting.
        return base;
    }
    return convertQuantityList(line->text(), combo->currentData().toDouble());
}
