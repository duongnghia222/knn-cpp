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
            throw out_of_range("Index out of range");
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
            throw out_of_range("Index out of range");
        }
        for (size_t i = index; i < size - 1; ++i) {
            array[i] = array[i + 1];
        }
        --size;
    }

    virtual T& get(int index) const override {
        if (index < 0 || index >= static_cast<int>(size)) {
            throw out_of_range("Index out of range");
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



class Dataset {
private:
    List<List<int>*>* data;
    List<string>* columnNames;

public:
    Dataset() : data(new DynamicArrayList<List<int>*>()), columnNames(new DynamicArrayList<string>()) {}

    ~Dataset() {
        for (int i = 0; i < data->length(); ++i) {
            delete data->get(i);
        }
        delete data;
        delete columnNames;
    }

    Dataset(const Dataset& other) : Dataset() {
        *this = other;
    }

    Dataset& operator=(const Dataset& other) {
        if (this != &other) {
            for (int i = 0; i < data->length(); ++i) {
                delete data->get(i);
            }
            data->clear();
            columnNames->clear();

            for (int i = 0; i < other.data->length(); ++i) {
                List<int>* newRow = new DynamicArrayList<int>();
                List<int>* otherRow = other.data->get(i);
                for (int j = 0; j < otherRow->length(); ++j) {
                    newRow->push_back(otherRow->get(j));
                }
                data->push_back(newRow);
            }

            for (int i = 0; i < other.columnNames->length(); ++i) {
                columnNames->push_back(other.columnNames->get(i));
            }
        }
        return *this;
    }

    bool loadFromCSV(const char* fileName) {
        ifstream file(fileName);
        if (!file.is_open()) {
            return false;
        }

        string line;
        // Read the column names
        if (getline(file, line)) {
            istringstream sstream(line);
            string columnName;
            while (getline(sstream, columnName, ',')) {
                columnNames->push_back(columnName);
            }
        } else {
            return false;
        }

        // Read the data
        while (getline(file, line)) {
            istringstream sstream(line);
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

    void printHead(int nRows = 5, int nCols = 5) const {
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

    void printTail(int nRows = 5, int nCols = 5) const {
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

    void getShape(int& nRows, int& nCols) const {
        nRows = data->length();
        nCols = data->length() > 0 ? data->get(0)->length() : 0;
    }

    void columns() const {
        for (int i = 0; i < columnNames->length(); ++i) {
            if (i > 0) cout << ", ";
            cout << columnNames->get(i);
        }
        cout << endl;
    }

    bool drop(int axis = 0, int index = 0, string columns = "") {
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

    Dataset extract(int startRow = 0, int endRow = -1, int startCol = 0, int endCol = -1) const {
        if (endRow == -1) {
            endRow = data->length();
        }
        if (endCol == -1) {
            endCol = data->get(0)->length();
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


    int length() const {
        return data->length();
    }

    // Returns the number of features (columns) in each item
    // Assuming all rows have the same number of columns
    int width() const {
        if (length() > 0) {
            return data->get(0)->length();
        }
        return 0;
    }

    // Returns the pointer to the data array at the given index (row)
    List<int>* get(int index) const {
        return data->get(index);
    }

    // Adds a new row to the dataset with the given data and size
    void add(List<int>* rowData) {
        List<int>* newRow = new DynamicArrayList<int>(); // Assuming DynamicArrayList is a concrete implementation of List
        for (int i = 0; i < rowData->length(); ++i) {
            newRow->push_back(rowData->get(i));
        }
        data->push_back(newRow);
}

    // Adds a new label to the dataset (for y_train and y_test)
    void add(int label) {
        List<int>* newRow = new DynamicArrayList<int>();
        newRow->push_back(label);
        data->push_back(newRow);
    }

    // Returns the maximum label value; used to create the count array for label frequencies
    int max_label() const {
        int maxLabel = -1;
        for (int i = 0; i < length(); ++i) {
            List<int>* row = get(i);
            if (row->length() > 0 && row->get(0) > maxLabel) {
                maxLabel = row->get(0);
            }
        }
        return maxLabel;
    }
};


double euclideanDistance(const List<int>* a, const List<int>* b, int size) {
    double distance = 0.0;
    for (int i = 0; i < size; ++i) {
        int diff = a->get(i) - b->get(i);
        distance += diff * diff;
    }
    return sqrt(distance);
}

class kNN {
private:
    int k;
    Dataset X_train;
    Dataset y_train;

public:
    kNN(int k = 5) : k(k) {}

    void fit(const Dataset& X_train, const Dataset& y_train) {
        this->X_train = X_train;
        this->y_train = y_train;
    }

    Dataset predict(const Dataset& X_test) {
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

    double score(const Dataset& y_test, const Dataset& y_pred) {
        int correctPredictions = 0;
        for (int i = 0; i < y_test.length(); ++i) {
            if (y_test.get(i) == y_pred.get(i)) {
                correctPredictions++;
            }
        }
        return static_cast<double>(correctPredictions) / y_test.length();
    }
};

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