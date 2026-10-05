#pragma once
#include "database_table_tab.h"

// Uses the same transactional Add/Modify editor as the other database tabs.
class ExperimentsTab : public DatabaseTableTab {
public:
    explicit ExperimentsTab(QWidget* parent = nullptr);
};
