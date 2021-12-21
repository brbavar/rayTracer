Color::Color() {}

Color::Color(const Color& orig) {
    std::cout << "begin Color copy" << std::endl;

    /* memcpy(&blue, &orig.blue, sizeof orig.blue);
    memcpy(&green, &orig.green, sizeof orig.green);
    memcpy(&red, &orig.red, sizeof orig.red);
    memcpy(&opacity, &orig.opacity, sizeof orig.opacity); */

    blue = orig.blue, green = orig.green, red = orig.red;
    opacity = orig.opacity;

    std::cout << "b = " << blue << std::endl;
    std::cout << "g = " << green << std::endl;
    std::cout << "r = " << red << std::endl;
    std::cout << "o = " << opacity << std::endl;

    std::cout << "end Color copy" << std::endl;
}

Color::Color(double b, double g, double r, double o) {
    std::cout << "begin Color" << std::endl;
    
    blue = b, green = g, red = r;
    opacity = o;

    std::cout << "b = " << blue << std::endl;
    std::cout << "g = " << green << std::endl;
    std::cout << "r = " << red << std::endl;
    std::cout << "o = " << opacity << std::endl;

    std::cout << "end Color" << std::endl;
}

Color::~Color() {
    std::cout << "destroy Color" << std::endl;
}

bool Color::operator==(Color other) {
    return blue == other.blue && green == other.green && red == other.red 
        && opacity == other.opacity;
}

bool Color::operator!=(Color other) {
    return !(*this == other);
}

Matrix::Matrix() {
    entries.reserve(1);
    entries[0].reserve(3);
}

Matrix::Matrix(const Matrix& orig) {
    std::cout << "begin Matrix copy" << std::endl;

    /* std::vector<double> list;
    entries.reserve(orig.entries.size());
    list.reserve(orig.entries[0].size());
    for(auto row : orig.entries) {
        for(double entry : row)
            list.push_back(entry);
        entries.push_back(list);
        list.clear();
    }
        
    memcpy(entries.data(), orig.entries.data(), sizeof orig.entries); */
    
    std::vector<double> list;
    entries.reserve(orig.entries.size());
    list.reserve(orig.entries[0].size()); 

    for(auto row : orig.entries) {
        for(double entry : row)
            list.push_back(entry);
        entries.push_back(list);
        list.clear();
    }

    for(auto row : entries)
        for(double entry : row)
            std::cout << entry << " ";

    std::cout << "\nend Matrix copy" << std::endl;
}

Matrix::Matrix(std::vector<std::vector<double> > items) {
    std::vector<double> list;
    entries.reserve(items.size());
    list.reserve(items[0].size());
    for(auto row : items) {
        for(double item : row)
            list.push_back(item);
        entries.push_back(list);
        list.clear();
    }
}

Matrix::Matrix(double x, double y) {
    entries.reserve(1);
    std::vector<double> coords;
    coords.reserve(2);
    coords.push_back(x);
    coords.push_back(y);
    entries.push_back(coords);
}

Matrix::Matrix(double x, double y, double z) {
    std::cout << "begin Matrix" << std::endl;

    entries.reserve(1);
    std::vector<double> coords;
    coords.reserve(3);
    coords.push_back(x);
    coords.push_back(y);
    coords.push_back(z);
    entries.push_back(coords);

    for(auto row : entries)
        for(double entry : row)
            std::cout << entry << " ";

    std::cout << "\nend Matrix" << std::endl;
}

Matrix::~Matrix() {
    std::cout << "destroy Matrix" << std::endl;
}

/* Returns -1 if either matrix has more than one row of entries, in which case one of them 
   cannot be said to represent a point in n-dimensional space. They must both be capable of 
   representing points in space, so that a distance between the points may be calculated. 
   Of course, there is no such thing as a negative distance. */
double Matrix::distanceTo(Matrix other) {
    if(entries.size() == 1 && other.entries.size() == 1) {
        double sumSqrs = 0;
        for(int c = 0; c < entries[0].size(); c++)
            sumSqrs += pow(other.entries[0][c] - entries[0][c], 2);
        return sqrt(sumSqrs);
    }
    return -1;
}

/* Returns -1 if matrix does not have exactly one row of entries, in which case it cannot
   be said to represent a vector (i.e., list of numbers rather than a table of numbers), 
   which has a magnitude. Of course, there is no such thing as a negative magnitude. */
double Matrix::mag() {
    double sumSqrs = 0;
    if(entries.size() == 1) {
        for(double e : entries[0])
            sumSqrs += e * e;
        return sqrt(sumSqrs);
    }
    return -1;
}

/* Returns empty matrix if this matrix has "magnitude" <= 0, because in that case
   this matrix has no magnitude with which we can normalize it. A nonzero magnitude
   is necessary for vector normalization, which is the kind of normalization we 
   are interested in. */
Matrix Matrix::normalize() {
    Matrix result = Matrix();
    double mag = this->mag();
    if(mag > 0) {
        std::vector<double> list;
        list.reserve(entries[0].size());
        for(auto row: entries) {
            for(double entry : row)
                list.push_back(entry / mag);
            result.entries.push_back(list);
            list.clear();
        }
    }
    return result;
}

bool Matrix::rowsColsEq(Matrix other) {
    bool match = entries.size() == other.entries.size();
    for(int r = 0; match && r < entries.size(); r++)
        match = entries[r].size() == other.entries[r].size();
    return match;
}

/* The dot product requires that the two matrices being multiplied have the same number of entries. 
   If they do not, then in the case of this function, a garbage value of 0 will be returned by default. */
double Matrix::dot(Matrix other) {
    double sumProds = 0;
    if(rowsColsEq(other))
        for(int r = 0; r < entries.size(); r++)
            for(int c = 0; c < entries[r].size(); c++)
                sumProds += entries[r][c] * other.entries[r][c];
    return sumProds;
}

/* For our purposes, and by convention, the cross product is only defined for coordinate vectors
   with exacty three entries. If this matrix and the other do not meet those criteria, then an
   empty matrix is returned. */
Matrix Matrix::cross(Matrix other) {
    Matrix result = Matrix();
    std::vector<double> list;
    if(entries.size() == 1 && other.entries.size() == 1 && 
        entries[0].size() == 3 && other.entries[0].size() == 3) {
        list.reserve(3);
        for(int i = 0; i < 3; i++)
            list.push_back(entries[0][i == 2 ? 0 : i + 1] * other.entries[0][i == 0 ? 2 : i - 1] 
                - entries[0][i == 0 ? 2 : i - 1] * other.entries[0][i == 2 ? 0 : i + 1]);
        result.entries.push_back(list);
    }
    return result;
}

Matrix Matrix::negate() {
    Matrix result = Matrix();
    std::vector<double> list;
    list.reserve(entries[0].size());
    for(auto row : entries) {
        for(double entry : row)
            list.push_back(-entry);
        result.entries.push_back(list);
        list.clear();
    }
    return result;
}

Matrix Matrix::operator*(double scalar) {
    Matrix result = Matrix();
    std::vector<double> list;
    list.reserve(entries[0].size());
    for(auto row : entries) {
        for(double entry : row)
            list.push_back(entry * scalar);
        result.entries.push_back(list);
        list.clear();
    }
    return result;
}

/* Call this matrix A and the other matrix B. In order to be multiplied, A must have exactly as
   many columns as B does rows. If this condition is unsatisfied, an empty matrix is returned. 
   Multiplying A and B gives us AB, which has exactly as many rows as A and exactly as many 
   columns as B. */
Matrix Matrix::operator*(Matrix other) {
    Matrix result = Matrix();
    if(entries[0].size() == other.entries.size()) {
        // Make room in AB for exactly as many rows as A has.
        result.entries.reserve(entries.size());

        /* Give every row in AB exactly as many entries as there are columns in B. 
        That way, AB will have exactly as many columns as B. */
        std::vector<double> row(other.entries[0].size(), 0);

        for(int r = 0; r < entries.size(); r++) {
            result.entries.push_back(row);
            for(int c = 0; c < result.entries[r].size(); c++)
                for(int n = 0; n < entries[0].size(); n++)
                    result.entries[r][c] += entries[r][n] * other.entries[n][c];
        }
    }
    return result;
}

/* Returns an empty matrix if the matrices being added do not meet the requirement of having 
   equal numbers of rows and columns */
Matrix Matrix::operator+(Matrix other) {
    Matrix result = Matrix();
    if(rowsColsEq(other)) {
        std::vector<double> list;
        list.reserve(entries[0].size());
        for(int r = 0; r < entries.size(); r++) {
            for(int c = 0; c < entries[r].size(); c++)
                list.push_back(entries[r][c] + other.entries[r][c]);
            result.entries.push_back(list);
            list.clear();
        }
    }
    return result;
}

Matrix Matrix::operator-(Matrix other) {
    return *this + other.negate();
}

bool Matrix::operator==(Matrix other) {
    if(rowsColsEq(other)) {
        for(int r = 0; r < entries.size(); r++)
            for(int c = 0; c < entries[r].size(); c++)
                if(entries[r][c] != other.entries[r][c])
                    return false;
        return true;
    }
    else 
        return false;
}

bool Matrix::operator!=(Matrix other) {
    return !(*this == other);
}

void Matrix::operator=(const Matrix& orig) {
    entries.clear();
    entries.reserve(orig.entries.size());
    std::vector<double> list;
    list.reserve(orig.entries[0].size());
    for(auto row : orig.entries) {
        for(double entry : row)
            list.push_back(entry);
        entries.push_back(list);
        list.clear();
    }
}


Light::Light() {}

Light::Light(double x, double y, double z, Color color) {
    std::cout << "begin Light" << std::endl;

    entries.reserve(1);
    std::vector<double> coords;                       // UNIT TEST THIS FUNC
    coords.reserve(3);
    coords.push_back(x);
    coords.push_back(y); 
    coords.push_back(z);
    entries.push_back(coords);

    clr = color;

    for(auto row : entries)
        for(double entry : row)
            std::cout << entry << " ";
    std::cout << "\nb = " << clr.blue << std::endl;
    std::cout << "g = " << clr.green << std::endl;
    std::cout << "r = " << clr.red << std::endl;
    std::cout << "o = " << clr.opacity << std::endl;

    std::cout << "end Light" << std::endl;
}

Ray::Ray() {}

/* Subtract coordinates of initial point from those of non-initial point,
    then normalize, and you get unit vector pointing same direction that ray points. */
Ray::Ray(Matrix initPt, Matrix nonInitPt) {
    entries.reserve(initPt.entries.size());
    std::vector<double> coords;
    coords.reserve(initPt.entries[0].size());
    for(auto row : initPt.entries) {
        for(double entry : row)
            coords.push_back(entry);
        entries.push_back(coords);
        coords.clear();
    }
    unitDir = (nonInitPt - initPt).normalize();
}

bool Ray::operator==(Ray other) {
    return *this == other && unitDir == other.unitDir;
}

bool Ray::operator!=(Ray other) {
    return !(*this == other);
}

Camera::Camera() {
    focus = Matrix(0, 2, 0);
}

Camera::Camera(Ray ray) {
    entries.reserve(ray.entries.size());
    std::vector<double> coords;
    coords.reserve(ray.entries[0].size());
    for(auto row : ray.entries) {
        for(double entry : row)
            coords.push_back(entry);            // UNIT TEST THIS FUNC
        entries.push_back(coords);
    }
    unitDir = ray.unitDir;
    right = Matrix(entries) + Matrix(1, 0, 0);
    focus = Matrix(0, 2, 0);
}

// Centered at origin of local coordinate system by default
Shape::Shape() {
    std::cout << "in Shape default" << std::endl;
    locCenter = Matrix(0, 0, 0);
}

Shape::Shape(const Shape& orig) {
    locCenter = orig.locCenter, globCenter = orig.globCenter;  // Some of these may have undefined vals,
    locRadius = orig.locRadius, globRadius = orig.globRadius;  // depending on the type of orig.
    locHeight = orig.locHeight, globHeight = orig.globHeight;  // MIGHT NEED TO CHANGE
    locWidth = orig.locWidth, globWidth = orig.globWidth;
    normal = orig.normal;
    type = orig.type;
    clr = orig.clr;
}

Shape::Shape(Color color) {
    clr = color;
}

Shape::Shape(Matrix c, double r, Color color) {
    std::cout << "in Shape" << std::endl;
    globCenter = c, globRadius = r;
    clr = color, type = "Sphere";
}

double Shape::intersects(Ray ray) {
    if(type == "Sphere") {
        // Solving for coefficients derived from equation for sphere
        double a = ray.unitDir.dot(ray.unitDir);
        Matrix initPt = ray - globCenter;
        double b = 2 * initPt.dot(ray.unitDir);
        double c = initPt.dot(initPt) - pow(globRadius, 2);

        if(b * b - 4 * a * c > 0) {
            double dist = (globCenter - ray).dot(ray.unitDir);
            Matrix ptBtwn = ray + ray.unitDir * dist;
            double offset = sqrt( pow(globRadius, 2) - pow((globCenter - ptBtwn).mag(), 2) );
            dist -= offset;
            return dist;
        }
    }
    if(type == "Plane") {
        double dist = ray.entries[0][1] / (-ray.unitDir.entries[0][1]);
        return dist;
    }
    return -1;
}

bool Shape::operator==(Shape other) {
    if(type != other.type || clr != other.clr)
        return false;

    bool sameCenter = true, sameDims = true, sameNormal = true;

    if(type == "Sphere" || type == "Plane") {
        sameCenter = /* locCenter == other.locCenter && */ globCenter == other.globCenter;
    }
    if(type == "Sphere") {
        sameDims = /* locRadius == other.locRadius && */ globRadius == other.globRadius;
    }
    if(type == "Plane") {
        sameDims = /* locHeight == other.locHeight && locWidth == other.locWidth 
            && */ globHeight == other.globHeight && globWidth == other.globWidth;
        sameNormal = normal == other.normal;
    }

    return sameCenter && sameDims && sameNormal;
}

bool Shape::operator!=(Shape other) {
    return !(*this == other);
}

void Shape::operator=(const Shape& orig) {
    locCenter = orig.locCenter, globCenter = orig.globCenter;  // Some of these may have undefined vals,
    locRadius = orig.locRadius, globRadius = orig.globRadius;  // depending on the type of orig.
    locHeight = orig.locHeight, globHeight = orig.globHeight;  // MIGHT NEED TO CHANGE
    locWidth = orig.locWidth, globWidth = orig.globWidth;
    normal = orig.normal;
    type = orig.type;
    clr = orig.clr;
}

Sphere::Sphere() {
    type = "Sphere";
}

Sphere::Sphere(Matrix c, double r, Color color) {
    std::cout << "in Sphere" << std::endl;
    globCenter = c, globRadius = r; 
    clr = color, type = "Sphere";
}

Plane::Plane() {
    type = "Plane";
}

Plane::Plane(Ray r, Matrix c, double h, double w) {
    normal = r, globCenter = c, globHeight = h, globWidth = w;
    double aspectRatio = w / h;
    locHeight = (w > h ? 1 / aspectRatio : 1);
    locWidth = (w > h ? 1 : aspectRatio);
    type = "Plane";
}

Plane::Plane(Ray r, Matrix c, double h, double w, Color color) {
    normal = r, globCenter = c, globHeight = h, globWidth = w;
    double aspectRatio = w / h;
    locHeight = (w > h ? 1 / aspectRatio : 1);
    locWidth = (w > h ? 1 : aspectRatio);
    clr = color, type = "Plane";
}