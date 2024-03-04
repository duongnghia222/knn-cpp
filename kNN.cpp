#include "kNN.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */

template<typename T>
class DynamicArrayList : public List<T> {
private:
    T* array;
    size_t capacity;
    size_t size;

    void expandCapacity() {
        capacity = capacity == 0 ? 1 : capacity * 2;
        T* newArray = new T[capacity];
        for (size_t i = 0; i < size; ++i) {
            newArray[i] = array[i];
        }
        delete[] array;
        array = newArray;
    }

public:
    DynamicArrayList() : array(nullptr), capacity(0), size(0) {}

    virtual ~DynamicArrayList() override {
        delete[] array;
    }

    virtual void push_back(T value) override {
        if (size == capacity) {
            expandCapacity();
        }
        array[size++] = value;
    }

    virtual void push_front(T value) override {
        if (size == capacity) {
            expandCapacity();
        }
        for (size_t i = size; i > 0; --i) {
            array[i] = array[i - 1];
        }
        array[0] = value;
        ++size;
    }

    virtual void insert(int index, T value) override {
        if (index < 0 || index > static_cast<int>(size)) {
            throw std::out_of_range("Index out of range");
        }
        if (size == capacity) {
            expandCapacity();
        }
        for (size_t i = size; i > static_cast<size_t>(index); --i) {
            array[i] = array[i - 1];
        }
        array[index] = value;
        ++size;
    }

    virtual void remove(int index) override {
        if (index < 0 || index >= static_cast<int>(size)) {
            throw std::out_of_range("Index out of range");
        }
        for (size_t i = index; i < size - 1; ++i) {
            array[i] = array[i + 1];
        }
        --size;
    }

    virtual T& get(int index) const override {
        if (index < 0 || index >= static_cast<int>(size)) {
            throw std::out_of_range("Index out of range");
        }
        return array[index];
    }

    virtual int length() const override {
        return static_cast<int>(size);
    }

    virtual void clear() override {
        delete[] array;
        array = nullptr;
        size = 0;
        capacity = 0;
    }

    virtual void print() const override {
        for (size_t i = 0; i < size; ++i) {
            cout << array[i] << ' ';
        }
        cout << endl;
    }

    virtual void reverse() override {
        for (size_t i = 0; i < size / 2; ++i) {
            T temp = array[i];
            array[i] = array[size - 1 - i];
            array[size - 1 - i] = temp;
        }
    }
};





    void Dataset::clearData() {
        if (data) {
            for (int i = 0; i < data->length(); ++i) {
                delete data->get(i);
            }
            delete data;
            data = nullptr; // Ensure the pointer is set to nullptr after deletion.
        }
    }

    Dataset::Dataset() : data(new DynamicArrayList<List<int>*>()), columnNames(new DynamicArrayList<string>()) {}

    Dataset::~Dataset() {
        clearData();
        delete columnNames;
    }

    Dataset::Dataset(const Dataset& other) : data(new DynamicArrayList<List<int>*>()), columnNames(new DynamicArrayList<string>()) {
        for (int i = 0; i < other.data->length(); ++i) {
            List<int>* newRow = new DynamicArrayList<int>();
            for (int j = 0; j < other.data->get(i)->length(); ++j) {
                newRow->push_back(other.data->get(i)->get(j));
            }
            data->push_back(newRow);
        }

        for (int i = 0; i < other.columnNames->length(); ++i) {
            columnNames->push_back(other.columnNames->get(i));
        }
    }

    // Assignment operator
    Dataset& Dataset::operator=(const Dataset& other) {
        if (this != &other) {
            // Create new data structures as a local copy.
            List<List<int>*>* newData = new DynamicArrayList<List<int>*>();
            List<string>* newColumnNames = new DynamicArrayList<string>();

            // Attempt to copy the data.
            for (int i = 0; i < other.data->length(); ++i) {
                List<int>* newRow = new DynamicArrayList<int>();
                for (int j = 0; j < other.data->get(i)->length(); ++j) {
                    newRow->push_back(other.data->get(i)->get(j));
                }
                newData->push_back(newRow);
            }
            for (int i = 0; i < other.columnNames->length(); ++i) {
                newColumnNames->push_back(other.columnNames->get(i));
            }

            // Delete the old data and assign the new data.
            clearData();
            delete columnNames;
            
            data = newData;
            columnNames = newColumnNames;
        }
        return *this;
    }

    bool Dataset::loadFromCSV(const char* fileName) {
        std::ifstream file(fileName);
        if (!file.is_open()) {
            return false;
        }

        string line;
        // Read the column names
        if (getline(file, line)) {
            std::istringstream sstream(line);
            string columnName;
            while (getline(sstream, columnName, ',')) {
                columnNames->push_back(columnName);
            }
        } else {
            return false;
        }

        // Read the data
        while (getline(file, line)) {
            std::istringstream sstream(line);
            string cell;
            List<int>* row = new DynamicArrayList<int>();
            
            while (getline(sstream, cell, ',')) {
                row->push_back(atoi(cell.c_str()));
            }
            data->push_back(row);
        }

        file.close();
        return true;
    }

    // Other methods inside the Dataset class:

    void Dataset::printHead(int nRows, int nCols) const {
        for (int col = 0; col < nCols && col < columnNames->length(); ++col) {
            if (col > 0) cout << ", ";
            cout << columnNames->get(col);
        }
        cout << endl;

        for (int row = 0; row < nRows && row < data->length(); ++row) {
            List<int>* rowData = data->get(row);
            for (int col = 0; col < nCols && col < rowData->length(); ++col) {
                if (col > 0) cout << ", ";
                cout << rowData->get(col);
            }
            cout << endl;
        }
    }

    void Dataset::printTail(int nRows, int nCols) const {
        int startRow = max(0, data->length() - nRows);
        for (int row = startRow; row < data->length(); ++row) {
            List<int>* rowData = data->get(row);
            for (int col = 0; col < nCols && col < rowData->length(); ++col) {
                if (col > 0) cout << ", ";
                cout << rowData->get(col);
            }
            cout << endl;
        }
    }

    void Dataset::getShape(int& nRows, int& nCols) const {
        nRows = data->length();
        nCols = data->length() > 0 ? data->get(0)->length() : 0;
    }

    void Dataset::columns() const {
        for (int i = 0; i < columnNames->length(); ++i) {
            if (i > 0) cout << ", ";
            cout << columnNames->get(i);
        }
        cout << endl;
    }

    bool Dataset::drop(int axis, int index, string columns) {
        if (axis == 0) {
            if (index < 0 || index >= data->length()) {
                return false;
            }
            delete data->get(index);
            data->remove(index);
        } else if (axis == 1) {
            int colIndex = -1;
            for (int i = 0; i < columnNames->length(); ++i) {
                if (columnNames->get(i) == columns) {
                    colIndex = i;
                    break;
                }
            }
            if (colIndex == -1) {
                return false;
            }
            for (int i = 0; i < data->length(); ++i) {
                List<int>* rowData = data->get(i);
                rowData->remove(colIndex);
            }
            columnNames->remove(colIndex);
        } else {
            return false;
        }
        return true;
    }

    Dataset Dataset::extract(int startRow, int endRow, int startCol, int endCol) const {
    // Adjust endRow and endCol to the length of the dataset if they are negative
        endRow = (endRow == -1) ? data->length() : endRow;
        endCol = (endCol == -1) ? (data->length() > 0 ? data->get(0)->length() : 0) : endCol;

        // Check if the specified ranges are within the valid range of the dataset dimensions
        if (startRow < 0 || startCol < 0 || endRow > data->length() || endCol > (data->length() > 0 ? data->get(0)->length() : 0)) {
            throw std::out_of_range("Specified range is out of the dataset dimensions");
        }

        Dataset extractedDataset;
        for (int col = startCol; col < endCol && col < columnNames->length(); ++col) {
            extractedDataset.columnNames->push_back(columnNames->get(col));
        }

        for (int row = startRow; row < endRow && row < data->length(); ++row) {
            List<int>* rowData = data->get(row);
            List<int>* extractedRow = new DynamicArrayList<int>();
            for (int col = startCol; col < endCol && col < rowData->length(); ++col) {
                extractedRow->push_back(rowData->get(col));
            }
            extractedDataset.data->push_back(extractedRow);
        }

        return extractedDataset;
}


    int Dataset::length() const {
        return data->length();
    }

    // Returns the number of features (columns) in each item
    // Assuming all rows have the same number of columns
    int Dataset::width() const {
        if (length() > 0) {
            return data->get(0)->length();
        }
        return 0;
    }

    // Returns the pointer to the data array at the given index (row)
    List<int>* Dataset::get(int index) const {
        return data->get(index);
    }

    // Adds a new row to the dataset with the given data and size
    void Dataset::add(List<int>* rowData) {
        List<int>* newRow = new DynamicArrayList<int>(); // Assuming DynamicArrayList is a concrete implementation of List
        for (int i = 0; i < rowData->length(); ++i) {
            newRow->push_back(rowData->get(i));
        }
        data->push_back(newRow);
}

    // Adds a new label to the dataset (for y_train and y_test)
    void Dataset::add(int label) {
        List<int>* newRow = new DynamicArrayList<int>();
        newRow->push_back(label);
        data->push_back(newRow);
    }

    // Returns the maximum label value; used to create the count array for label frequencies
    int Dataset::max_label() const {
        int maxLabel = -1;
        for (int i = 0; i < length(); ++i) {
            List<int>* row = get(i);
            if (row->length() > 0 && row->get(0) > maxLabel) {
                maxLabel = row->get(0);
            }
        }
        return maxLabel;
    }



double euclideanDistance(const List<int>* a, const List<int>* b, int size) {
    double distance = 0.0;
    for (int i = 0; i < size; ++i) {
        int diff = a->get(i) - b->get(i);
        distance += diff * diff;
    }
    return std::sqrt(distance);
}

    kNN::kNN(int k) : k(k) {}

    void kNN::fit(const Dataset& X_train, const Dataset& y_train) {
        this->X_train = X_train;
        this->y_train = y_train;
    }

    Dataset kNN::predict(const Dataset& X_test) {
        Dataset predictions;

        for (int i = 0; i < X_test.length(); ++i) {
            double* distances = new double[X_train.length()];
            int* indices = new int[X_train.length()];

            for (int j = 0; j < X_train.length(); ++j) {
                distances[j] = euclideanDistance(X_test.get(i), X_train.get(j), X_train.width());
                indices[j] = j;
            }

            // Simple selection sort to find the k smallest distances
            for (int j = 0; j < k; ++j) {
                int minIndex = j;
                for (int l = j + 1; l < X_train.length(); ++l) {
                    if (distances[l] < distances[minIndex]) {
                        minIndex = l;
                    }
                }
                swap(distances[j], distances[minIndex]);
                swap(indices[j], indices[minIndex]);
            }

            // Majority voting for the nearest k neighbors
            int* count = new int[y_train.max_label() + 1](); // Assuming labels are 0-indexed
            for (int j = 0; j < k; ++j) {
                int label = y_train.get(indices[j])->get(0);
                count[label]++;
            }

            int maxCount = 0;
            int predictedLabel = -1;
            for (int j = 0; j <= y_train.max_label(); ++j) {
                if (count[j] > maxCount) {
                    maxCount = count[j];
                    predictedLabel = j;
                }
            }

            predictions.add(predictedLabel);

            delete[] distances;
            delete[] indices;
            delete[] count;
        }

        return predictions;
    }

    double kNN::score(const Dataset& y_test, const Dataset& y_pred) {
    int correctPredictions = 0;
    for (int i = 0; i < y_test.length(); ++i) {
        // Compare the values instead of pointers
        if (y_test.get(i)->get(0) == y_pred.get(i)->get(0)) {
            correctPredictions++;
        }
    }
    return static_cast<double>(correctPredictions) / y_test.length();
}


void train_test_split(Dataset& X, Dataset& y, double test_size, 
                      Dataset& X_train, Dataset& X_test, Dataset& y_train, Dataset& y_test) {
    int totalRows = X.length();
    int trainRows = static_cast<int>((1 - test_size) * totalRows);
    int testRows = totalRows - trainRows;

    for (int i = 0; i < trainRows; ++i) {
        X_train.add(X.get(i));
        y_train.add(y.get(i));
    }
    for (int i = trainRows; i < totalRows; ++i) {
        X_test.add(X.get(i));
        y_test.add(y.get(i));
    }
}