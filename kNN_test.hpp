#include "main.hpp"

/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */

template <typename T>
class List
{
public:
    virtual ~List() = default;
    virtual void push_back(T value) = 0;
    virtual void push_front(T value) = 0;
    virtual void insert(int index, T value) = 0;
    virtual void remove(int index) = 0;
    virtual T &get(int index) const = 0;
    virtual int length() const = 0;
    virtual void clear() = 0;
    virtual void print() const = 0;
    virtual void reverse() = 0;
};

template <typename T>
struct Node
{
    T data;
    Node *next;

    Node(const T &value) : data(value), next(nullptr) {}
};

template <typename T>
class DynamicLinkedList : public List<T>
{
private:
    Node<T> *head;
    size_t size;

public:
    DynamicLinkedList() : head(nullptr), size(0) {}

    ~DynamicLinkedList() override
    {
        clear();
    }

    void push_back(T value) override
    {
        if (!head)
        {
            head = new Node<T>(value);
        }
        else
        {
            Node<T> *current = head;
            while (current->next)
            {
                current = current->next;
            }
            current->next = new Node<T>(value);
        }
        ++size;
    }

    void push_front(T value) override
    {
        Node<T> *newNode = new Node<T>(value);
        newNode->next = head;
        head = newNode;
        ++size;
    }

    void insert(int index, T value) override
    {
        if (index < 0 || index > static_cast<int>(size))
        {
            cout << "Index out of range";
            return;
        }

        if (index == 0)
        {
            push_front(value);
        }
        else if (index == size)
        {
            push_back(value);
        }
        else
        {
            Node<T> *current = head;
            for (int i = 0; i < index - 1; ++i)
            {
                current = current->next;
            }
            Node<T> *newNode = new Node<T>(value);
            newNode->next = current->next;
            current->next = newNode;
            ++size;
        }
    }

    void remove(int index) override
    {
        if (index < 0 || index >= static_cast<int>(size))
        {
            cout << "Index out of range";
            return;
        }

        Node<T> *temp;
        if (index == 0)
        {
            temp = head;
            head = head->next;
        }
        else
        {
            Node<T> *current = head;
            for (int i = 0; i < index - 1; ++i)
            {
                current = current->next;
            }
            temp = current->next;
            current->next = temp->next;
        }
        delete temp;
        --size;
    }

    T &get(int index) const override
    {
        if (index < 0 || index >= static_cast<int>(size))
        {
            cout << "Index out of range";
            // Returning a reference to a temporary is bad practice,
            // but this is just a placeholder.
            return head->data;
        }

        Node<T> *current = head;
        for (int i = 0; i < index; ++i)
        {
            current = current->next;
        }
        return current->data;
    }

    int length() const override
    {
        return static_cast<int>(size);
    }

    void clear() override
    {
        Node<T> *current = head;
        while (current)
        {
            Node<T> *next = current->next;
            delete current;
            current = next;
        }
        head = nullptr;
        size = 0;
    }

    void print() const override
    {
        Node<T> *current = head;
        while (current)
        {
            cout << current->data << " ";
            current = current->next;
        }
        cout << endl;
    }

    void reverse() override
    {
        if (!head || !head->next)
        {
            return;
        }

        Node<T> *prev = nullptr;
        Node<T> *current = head;
        Node<T> *nextNode;

        while (current)
        {
            nextNode = current->next;
            current->next = prev;
            prev = current;
            current = nextNode;
        }

        head = prev;
    }
};

class Dataset
{
private:
    List<List<int> *> *data;
    List<string> *columnNames;
    void clearData()
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
    // You may need to define more
public:
    Dataset() : data(new DynamicLinkedList<List<int> *>()), columnNames(new DynamicLinkedList<string>()) {}
    ~Dataset()
    {
        clearData();
        delete columnNames;
    }
    Dataset(const Dataset &other) : data(new DynamicLinkedList<List<int> *>()), columnNames(new DynamicLinkedList<string>())
    {
        for (int i = 0; i < other.data->length(); ++i)
        {
            List<int> *newRow = new DynamicLinkedList<int>();
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
    Dataset &operator=(const Dataset &other)
    {
        if (this != &other)
        {
            // Create new data structures as a local copy.
            List<List<int> *> *newData = new DynamicLinkedList<List<int> *>();
            List<string> *newColumnNames = new DynamicLinkedList<string>();

            // Attempt to copy the data.
            for (int i = 0; i < other.data->length(); ++i)
            {
                List<int> *newRow = new DynamicLinkedList<int>();
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
    bool loadFromCSV(const char *fileName)
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
            List<int> *row = new DynamicLinkedList<int>();

            while (getline(sstream, cell, ','))
            {
                row->push_back(atoi(cell.c_str()));
            }
            data->push_back(row);
        }

        file.close();
        return true;
    }
    void printHead(int nRows = 5, int nCols = 5) const
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
        if (data->length() > 0)
        {
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
            if (row < nRows - 1 && row < data->length() - 1)
            {
                cout << endl;
            }
        }
    }
    void printTail(int nRows = 5, int nCols = 5) const
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
        if (totalRows > 0)
        {
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
            if (row < totalRows - 1)
            {
                cout << endl;
            }
        }
    }
    void getShape(int &nRows, int &nCols) const
    {
        nRows = data->length();
        nCols = data->length() > 0 ? data->get(0)->length() : 0;
    }
    void columns() const
    {
        for (int i = 0; i < columnNames->length(); ++i)
        {
            if (i > 0)
                cout << ", ";
            cout << columnNames->get(i);
        }
        cout << endl;
    }
    bool drop(int axis = 0, int index = 0, std::string columns = "")
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
    Dataset extract(int startRow = 0, int endRow = -1, int startCol = 0, int endCol = -1) const
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
            List<int> *extractedRow = new DynamicLinkedList<int>();
            for (int col = startCol; col <= endCol && col < rowData->length(); ++col)
            {
                extractedRow->push_back(rowData->get(col));
            }
            extractedDataset.data->push_back(extractedRow);
        }

        return extractedDataset;
    }
    int length() const
    {
        return data->length();
    }
    int width() const
    {
        if (length() > 0)
        {
            return data->get(0)->length();
        }
        return 0;
    }
    List<int> *get(int index) const
    {
        List<int> *row = data->get(index);
        return row;
    }
    void add(List<int> *rowData)
    {
        List<int> *newRow = new DynamicLinkedList<int>();
        for (int i = 0; i < rowData->length(); ++i)
        {
            newRow->push_back(rowData->get(i));
        }
        data->push_back(newRow);
    }
    void add(int label)
    {
        List<int> *newRow = new DynamicLinkedList<int>();
        newRow->push_back(label);
        data->push_back(newRow);
    }
    int max_label() const
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
    List<List<int> *> *getData() const
    {
        return data;
    }
    void addColumnNames(const List<string> &names)
    {
        columnNames->clear();
        for (int i = 0; i < names.length(); ++i)
        {
            columnNames->push_back(names.get(i));
        }
    }
    List<string> *getColumnNames() const
    {
        return columnNames;
    }
};

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

void quickSort(double *distances, int *indices, int left, int right)
{
    if (left >= right)
        return;

    double pivot = distances[(left + right) / 2];
    int i = left;
    int j = right;

    while (i <= j)
    {
        while (distances[i] < pivot)
            i++;
        while (distances[j] > pivot)
            j--;
        if (i <= j)
        {
            std::swap(distances[i], distances[j]);
            std::swap(indices[i], indices[j]);
            i++;
            j--;
        }
    }

    quickSort(distances, indices, left, j);
    quickSort(distances, indices, i, right);
}

class kNN
{
private:
    int k;
    Dataset X_train;
    Dataset y_train;

    // You may need to define more
public:
    kNN(int k = 5)
    {
        this->k = k;
    }
    void fit(const Dataset &X_train, const Dataset &y_train)
    {
        this->X_train = X_train;
        this->y_train = y_train;
    }
    Dataset predict(const Dataset &X_test)
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
                List<int> *row = y_train.get(indices[j]);
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
    double score(const Dataset &y_test, const Dataset &y_pred)
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
};

void train_test_split(Dataset &X, Dataset &y, double test_size,
                      Dataset &X_train, Dataset &X_test, Dataset &y_train, Dataset &y_test)
{
    int totalRows = X.length();
    int testRows = static_cast<int>(ceil(test_size * totalRows));
    int trainRows = totalRows - testRows;

    List<string> *X_columnNames = X.getColumnNames();
    List<string> *y_columnNames = y.getColumnNames();

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

double euclideanDistance(const List<int> *a, const List<int> *b, int size);
// Please add more or modify as needed
