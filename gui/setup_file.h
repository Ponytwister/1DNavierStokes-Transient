#pragma once
#include "model_controls.h"

namespace setup_file {
struct Setup {
    QString database, outputDirectory;
    std::vector<model_controls::Row> controls;
};
// Paths are absolute when saved; relative paths in files resolve beside the file.
void save(const QString& filename, const Setup& setup);
Setup load(const QString& filename);
// Restores values transactionally, requiring the same control names in the database.
void restore(const Setup& setup, const model_controls::Snapshot& original);
}
