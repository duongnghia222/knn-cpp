#include "kNN.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */


template <typename T>
struct Node {
    T data;
    Node* next;

    Node(const T& value) : data(value), next(nullptr) {}
};

template <typename T>
class DynamicLinkedList : public List<T> {
private:
    Node<T>* head;
    size_t size;

public:
    DynamicLinkedList() : head(nullptr), size(0) {}

    ~DynamicLinkedList() override {
        clear();
    }

    void push_back(const T& value) override {
        if (!head) {
            head = new Node<T>(value);
        } else {
            Node<T>* current = head;
            while (current->next) {
                current = current->next;
            }
            current->next = new Node<T>(value);
        }
        ++size;
    }

    void push_front(const T& value) override {
        Node<T>* newNode = new Node<T>(value);
        newNode->next = head;
        head = newNode;
        ++size;
    }

    void insert(int index, const T& value) override {
        if (index < 0 || index > static_cast<int>(size)) {
            cout << "Index out of range";
            return;
        }

        if (index == 0) {
            push_front(value);
        } else if (index == size) {
            push_back(value);
        } else {
            Node<T>* current = head;
            for (int i = 0; i < index - 1; ++i) {
                current = current->next;
            }
            Node<T>* newNode = new Node<T>(value);
            newNode->next = current->next;
            current->next = newNode;
            ++size;
        }
    }

    void remove(int index) override {
        if (index < 0 || index >= static_cast<int>(size)) {
            cout << "Index out of range";
            return;
        }

        Node<T>* temp;
        if (index == 0) {
            temp = head;
            head = head->next;
        } else {
            Node<T>* current = head;
            for (int i = 0; i < index - 1; ++i) {
                current = current->next;
            }
            temp = current->next;
            current->next = temp->next;
        }
        delete temp;
        --size;
    }

    T& get(int index) const override {
        if (index < 0 || index >= static_cast<int>(size)) {
            cout << "Index out of range";
            // Returning a reference to a temporary is bad practice,
            // but this is just a placeholder.
            return head->data;
        }

        Node<T>* current = head;
        for (int i = 0; i < index; ++i) {
            current = current->next;
        }
        return current->data;
    }

    int length() const override {
        return static_cast<int>(size);
    }

    void clear() override {
        Node<T>* current = head;
        while (current) {
            Node<T>* next = current->next;
            delete current;
            current = next;
        }
        head = nullptr;
        size = 0;
    }

    void print() const override {
        Node<T>* current = head;
        while (current) {
            cout << current->data << " ";
            current = current->next;
        }
        cout << endl;
    }

    void reverse() override {
        if (!head || !head->next) {
            return;
        }

        Node<T>* prev = nullptr;
        Node<T>* current = head;
        Node<T>* nextNode;

        while (current) {
            nextNode = current->next;
            current->next = prev;
            prev = current;
            current = nextNode;
        }

        head = prev;
    }
};




template <typename T>
class DynamicArrayList : public List<T>
{
private:
    T *array;
    size_t capacity;
    size_t size;

    void expandCapacity()
    {
        capacity = capacity == 0 ? 1 : capacity * 2;
        T *newArray = new T[capacity];
        for (size_t i = 0; i < size; ++i)
        {
            newArray[i] = array[i];
        }
        delete[] array;
        array = newArray;
    }

public:
    DynamicArrayList() : array(nullptr), capacity(0), size(0) {}

    virtual ~DynamicArrayList() override
    {
        delete[] array;
    }

    virtual void push_back(T value) override
    {
        if (size == capacity)
        {
            expandCapacity();
        }
        array[size++] = value;
    }

    virtual void push_front(T value) override
    {
        if (size == capacity)
        {
            expandCapacity();
        }
        for (size_t i = size; i > 0; --i)
        {
            array[i] = array[i - 1];
        }
        array[0] = value;
        ++size;
    }

    virtual void insert(int index, T value) override
    {
        if (index < 0 || index > static_cast<int>(size))
        {
            cout << ("Index out of range");
        }
        if (size == capacity)
        {
            expandCapacity();
        }
        for (size_t i = size; i > static_cast<size_t>(index); --i)
        {
            array[i] = array[i - 1];
        }
        array[index] = value;
        ++size;
    }

    virtual void remove(int index) override
    {
        if (index < 0 || index >= static_cast<int>(size))
        {
            cout << ("Index out of range");
        }
        for (size_t i = index; i < size - 1; ++i)
        {
            array[i] = array[i + 1];
        }
        --size;
    }

    virtual T &get(int index) const override
    {
        if (index < 0 || index >= static_cast<int>(size))
        {
            cout << ("Index out of range");
        }
        return array[index];
    }

    virtual int length() const override
    {
        return static_cast<int>(size);
    }

    virtual void clear() override
    {
        delete[] array;
        array = nullptr;
        size = 0;
        capacity = 0;
    }

    virtual void print() const override
    {
        for (size_t i = 0; i < size; ++i)
        {
            cout << array[i] << ' ';
        }
        cout << endl;
    }

    virtual void reverse() override
    {
        for (size_t i = 0; i < size / 2; ++i)
        {
            T temp = array[i];
            array[i] = array[size - 1 - i];
            array[size - 1 - i] = temp;
        }
    }
};

void Dataset::clearData()
{
    if (data)
    {
        for (int i = 0; i < data->length(); ++i)
        {
            delete data->get(i);
        }
        delete data;
        data = nullptr; 
    }
}

Dataset::Dataset() : data(new DynamicArrayList<List<int> *>()), columnNames(new DynamicArrayList<string>()) {}

Dataset::~Dataset()
{
    clearData();
    delete columnNames;
}

Dataset::Dataset(const Dataset &other) : data(new DynamicArrayList<List<int> *>()), columnNames(new DynamicArrayList<string>())
{
    for (int i = 0; i < other.data->length(); ++i)
    {
        List<int> *newRow = new DynamicArrayList<int>();
        for (int j = 0; j < other.data->get(i)->length(); ++j)
        {
            newRow->push_back(other.data->get(i)->get(j));
        }
        data->push_back(newRow);
    }

    for (int i = 0; i < other.columnNames->length(); ++i)
    {
        columnNames->push_back(other.columnNames->get(i));
    }
}

// Assignment operator
Dataset &Dataset::operator=(const Dataset &other)
{
    if (this != &other)
    {
        // Create new data structures as a local copy.
        List<List<int> *> *newData = new DynamicArrayList<List<int> *>();
        List<string> *newColumnNames = new DynamicArrayList<string>();

        // Attempt to copy the data.
        for (int i = 0; i < other.data->length(); ++i)
        {
            List<int> *newRow = new DynamicArrayList<int>();
            for (int j = 0; j < other.data->get(i)->length(); ++j)
            {
                newRow->push_back(other.data->get(i)->get(j));
            }
            newData->push_back(newRow);
        }
        for (int i = 0; i < other.columnNames->length(); ++i)
        {
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

bool Dataset::loadFromCSV(const char *fileName)
{
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        return false;
    }

    string line;
    if (getline(file, line))
    {
        std::istringstream sstream(line);
        string columnName;
        while (getline(sstream, columnName, ','))
        {
            columnNames->push_back(columnName);
        }
    }
    else
    {
        return false;
    }

    while (getline(file, line))
    {
        std::istringstream sstream(line);
        string cell;
        List<int> *row = new DynamicArrayList<int>();

        while (getline(sstream, cell, ','))
        {
            row->push_back(atoi(cell.c_str()));
        }
        data->push_back(row);
    }

    file.close();
    return true;
}


void Dataset::printHead(int nRows, int nCols) const
{
    if (nRows < 0 || nCols < 0)
        return;

    // Print column names
    for (int col = 0; col < nCols && col < columnNames->length(); ++col)
    {
        if (col > 0)
            cout << " ";
        cout << columnNames->get(col);
    }
    if (data->length() > 0) {
        cout << endl;
    }

    // Print data rows
    for (int row = 0; row < nRows && row < data->length(); ++row)
    {
        List<int> *rowData = data->get(row);
        for (int col = 0; col < nCols && col < rowData->length(); ++col)
        {
            if (col > 0)
                cout << " ";
            cout << rowData->get(col);
        }
        if (row < nRows - 1 && row < data->length() - 1) {
            cout << endl;
        }
    }
}

void Dataset::printTail(int nRows, int nCols) const
{
    if (nRows < 0 || nCols < 0)
        return;

    int totalRows = data->length();
    int totalCols = columnNames->length();

    nRows = min(nRows, totalRows);
    nCols = min(nCols, totalCols);

    int startRow = max(0, totalRows - nRows);

    for (int col = max(0, totalCols - nCols); col < totalCols; ++col)
    {
        cout << columnNames->get(col);
        if (col < totalCols - 1)
            cout << " ";
    }
    if (totalRows > 0) {
        cout << endl;
    }

    for (int row = startRow; row < totalRows; ++row)
    {
        List<int> *rowData = data->get(row);
        for (int col = max(0, rowData->length() - nCols); col < rowData->length(); ++col)
        {
            cout << rowData->get(col);
            if (col < rowData->length() - 1)
                cout << " ";
        }
        if (row < totalRows - 1) {
            cout << endl;
        }
    }
}

void Dataset::getShape(int &nRows, int &nCols) const
{
    nRows = data->length();
    nCols = data->length() > 0 ? data->get(0)->length() : 0;
}

void Dataset::columns() const
{
    for (int i = 0; i < columnNames->length(); ++i)
    {
        if (i > 0)
            cout << ", ";
        cout << columnNames->get(i);
    }
    cout << endl;
}

bool Dataset::drop(int axis, int index, string columns)
{
    if (axis == 0)
    {
        if (index < 0 || index >= data->length())
        {
            return false;
        }
        delete data->get(index);
        data->remove(index);
    }
    else if (axis == 1)
    {
        int colIndex = -1;
        for (int i = 0; i < columnNames->length(); ++i)
        {
            if (columnNames->get(i) == columns)
            {
                colIndex = i;
                break;
            }
        }
        if (colIndex == -1)
        {
            return false;
        }
        for (int i = 0; i < data->length(); ++i)
        {
            List<int> *rowData = data->get(i);
            rowData->remove(colIndex);
        }
        columnNames->remove(colIndex);
    }
    else
    {
        return false;
    }
    return true;
}

Dataset Dataset::extract(int startRow, int endRow, int startCol, int endCol) const
{
    endRow = (endRow == -1) ? data->length() - 1 : endRow;
    endCol = (endCol == -1) ? (data->length() > 0 ? data->get(0)->length() - 1 : 0) : endCol;

    if (startRow < 0 || startCol < 0 || endRow >= data->length() || endCol >= (data->length() > 0 ? data->get(0)->length() : 0))
    {
        cout << "Specified range is out of the dataset dimensions" << endl;
        return Dataset();
    }

    Dataset extractedDataset;
    for (int col = startCol; col <= endCol && col < columnNames->length(); ++col)
    {
        extractedDataset.columnNames->push_back(columnNames->get(col));
    }

    for (int row = startRow; row <= endRow && row < data->length(); ++row)
    {
        List<int> *rowData = data->get(row);
        List<int> *extractedRow = new DynamicArrayList<int>();
        for (int col = startCol; col <= endCol && col < rowData->length(); ++col)
        {
            extractedRow->push_back(rowData->get(col));
        }
        extractedDataset.data->push_back(extractedRow);
    }

    return extractedDataset;
}


int Dataset::length() const
{
    return data->length();
}

int Dataset::width() const
{
    if (length() > 0)
    {
        return data->get(0)->length();
    }
    return 0;
}

List<int> *Dataset::get(int index) const
{
    List<int> * row = data->get(index);
    return row;
}

void Dataset::add(List<int> *rowData)
{
    List<int> *newRow = new DynamicArrayList<int>(); 
    for (int i = 0; i < rowData->length(); ++i)
    {
        newRow->push_back(rowData->get(i));
    }
    data->push_back(newRow);
}

void Dataset::add(int label)
{
    List<int> *newRow = new DynamicArrayList<int>();
    newRow->push_back(label);
    data->push_back(newRow);
}

int Dataset::max_label() const
{
    int maxLabel = -1;
    for (int i = 0; i < length(); ++i)
    {
        List<int> *row = get(i);
        if (row->length() > 0 && row->get(0) > maxLabel)
        {
            maxLabel = row->get(0);
        }
    }
    return maxLabel;
}

List<List<int>*>* Dataset::getData() const {
    return data;
}

void Dataset::addColumnNames(const List<string>& names) {
    columnNames->clear();
    for (int i = 0; i < names.length(); ++i) {
        columnNames->push_back(names.get(i));
    }
}

List<string>* Dataset::getColumnNames() const {
    return columnNames;
}


double euclideanDistance(const List<int> *a, const List<int> *b, int size)
{
    double distance = 0.0;
    for (int i = 0; i < size; ++i)
    {
        int diff = a->get(i) - b->get(i);
        distance += diff * diff;
    }
    return std::sqrt(distance);
}


void quickSort(double *distances, int *indices, int left, int right) {
    if (left >= right) return;
    
    double pivot = distances[(left + right) / 2];
    int i = left;
    int j = right;
    
    while (i <= j) {
        while (distances[i] < pivot) i++;
        while (distances[j] > pivot) j--;
        if (i <= j) {
            std::swap(distances[i], distances[j]);
            std::swap(indices[i], indices[j]);
            i++;
            j--;
        }
    }
    
    quickSort(distances, indices, left, j);
    quickSort(distances, indices, i, right);
}


kNN::kNN(int k) : k(k) {}

void kNN::fit(const Dataset &X_train, const Dataset &y_train)
{
    this->X_train = X_train;
    this->y_train = y_train;
}

Dataset kNN::predict(const Dataset &X_test)
{
    Dataset predictions;
    predictions.addColumnNames(*y_train.getColumnNames());
    for (int i = 0; i < X_test.length(); ++i)
    {
        double *distances = new double[X_train.length()];
        int *indices = new int[X_train.length()];

        for (int j = 0; j < X_train.length(); ++j)
        {
            distances[j] = euclideanDistance(X_test.get(i), X_train.get(j), X_train.width());
            indices[j] = j;
        }

        quickSort(distances, indices, 0, X_train.length() - 1);


        // Majority voting for the nearest k neighbors
        int *count = new int[y_train.max_label() + 1](); 
        for (int j = 0; j < k; ++j)
        {
            List<int> * row = y_train.get(indices[j]);
            int label = row->get(0);
            count[label]++;
        }

        int maxCount = 0;
        int predictedLabel = -1;
        for (int j = 0; j <= y_train.max_label(); ++j)
        {
            if (count[j] > maxCount)
            {
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

double kNN::score(const Dataset &y_test, const Dataset &y_pred)
{
    int correctPredictions = 0;
    for (int i = 0; i < y_test.length(); ++i)
    {
        // Compare the values instead of pointers
        if (y_test.get(i)->get(0) == y_pred.get(i)->get(0))
        {
            correctPredictions++;
        }
    }
    return static_cast<double>(correctPredictions) / y_test.length();
}

void train_test_split(Dataset &X, Dataset &y, double test_size,
                      Dataset &X_train, Dataset &X_test, Dataset &y_train, Dataset &y_test)
{
    int totalRows = X.length();
    int testRows = static_cast<int>(ceil(test_size * totalRows));
    int trainRows = totalRows - testRows;

    List<string>* X_columnNames = X.getColumnNames();
    List<string>* y_columnNames = y.getColumnNames();

    // Add column names to X_train and X_test
    X_train.addColumnNames(*X_columnNames);
    X_test.addColumnNames(*X_columnNames);

    // Add column names to y_train and y_test
    y_train.addColumnNames(*y_columnNames);
    y_test.addColumnNames(*y_columnNames);

    for (int i = 0; i < trainRows; ++i)
    {
        X_train.add(X.get(i));
        y_train.add(y.get(i));
    }


    for (int i = trainRows; i < totalRows; ++i)
    {
        X_test.add(X.get(i));
        y_test.add(y.get(i));
    }
}
