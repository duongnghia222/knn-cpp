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
            std::cout << array[i] << ' ';
        }
        std::cout << std::endl;
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
    List<std::string>* columnNames;

public:
    Dataset() : data(new DynamicArrayList<List<int>*>()), columnNames(new DynamicArrayList<std::string>()) {}

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
        std::ifstream file(fileName);
        if (!file.is_open()) {
            return false;
        }

        std::string line;
        // Read the column names
        if (getline(file, line)) {
            std::istringstream sstream(line);
            std::string columnName;
            while (getline(sstream, columnName, ',')) {
                columnNames->push_back(columnName);
            }
        } else {
            return false;
        }

        // Read the data
        while (getline(file, line)) {
            std::istringstream sstream(line);
            std::string cell;
            List<int>* row = new DynamicArrayList<int>();
            
            while (getline(sstream, cell, ',')) {
                row->push_back(std::atoi(cell.c_str()));
            }
            data->push_back(row);
        }

        file.close();
        return true;
    }

    // Other methods inside the Dataset class:

    void printHead(int nRows = 5, int nCols = 5) const {
        for (int col = 0; col < nCols && col < columnNames->length(); ++col) {
            if (col > 0) std::cout << ", ";
            std::cout << columnNames->get(col);
        }
        std::cout << std::endl;

        for (int row = 0; row < nRows && row < data->length(); ++row) {
            List<int>* rowData = data->get(row);
            for (int col = 0; col < nCols && col < rowData->length(); ++col) {
                if (col > 0) std::cout << ", ";
                std::cout << rowData->get(col);
            }
            std::cout << std::endl;
        }
    }

    void printTail(int nRows = 5, int nCols = 5) const {
        int startRow = std::max(0, data->length() - nRows);
        for (int row = startRow; row < data->length(); ++row) {
            List<int>* rowData = data->get(row);
            for (int col = 0; col < nCols && col < rowData->length(); ++col) {
                if (col > 0) std::cout << ", ";
                std::cout << rowData->get(col);
            }
            std::cout << std::endl;
        }
    }

    void getShape(int& nRows, int& nCols) const {
        nRows = data->length();
        nCols = data->length() > 0 ? data->get(0)->length() : 0;
    }

    void columns() const {
        for (int i = 0; i < columnNames->length(); ++i) {
            if (i > 0) std::cout << ", ";
            std::cout << columnNames->get(i);
        }
        std::cout << std::endl;
    }

    bool drop(int axis = 0, int index = 0, std::string columns = "") {
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
};


