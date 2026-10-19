#pragma once
#include <QWidget>
#include <QString>
#include <functional>
#include <optional>
#include <vector>
#include <string>

class QTableWidget;
class QPlainTextEdit;
class QPushButton;

class RunHistoryTab : public QWidget {
public:
    explicit RunHistoryTab(QWidget* parent = nullptr);
    void setDatabase(const QString& database);
    std::vector<std::pair<std::string, std::optional<std::string>>> selectedControls() const;
    std::function<void()> loadRequested;
private:
    void refresh();
    QString database_;
    QTableWidget* runs_;
    QPlainTextEdit* details_;
    QPushButton* load_;
    qint64 selectedRun_ = 0;
    std::vector<std::pair<std::string, std::optional<std::string>>> controls_;
};
