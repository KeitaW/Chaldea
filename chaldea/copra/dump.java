Compiled from "COPRA.java"
public class COPRA {
  static java.lang.String separator;
  static java.lang.String welcome1;
  static java.lang.String welcome2;
  static int delVerts;
  static int delClusters;
  static int simpleDelVerts;
  static int simpleDelClusters;
  public COPRA();
  public static void main(java.lang.String[]);
  public static VecPair clusterGraph(int, java.lang.String, java.lang.String, java.lang.String, java.lang.String, float, int, boolean, boolean, java.lang.String, int, boolean, boolean);
  public static VecPair clusterGraph(int, java.lang.String, java.lang.String, java.lang.String, java.lang.String, float, int, boolean, boolean, java.lang.String, int, boolean);
  public static void propagate(java.lang.String, int, java.lang.String, java.lang.String, java.lang.String, float);
  public static void readBiGraphEdges(java.lang.String, java.util.List<java.util.HashSet<java.lang.Integer>>, java.util.List<java.util.HashSet<java.lang.Integer>>, java.util.List<java.lang.String>, java.util.List<java.lang.String>, java.util.List<java.util.HashMap<java.lang.Integer, java.lang.Float>>);
  static {};
}
Compiled from "COPRA.java"
class ClusterLabel {
  public ClusterLabel(int, float);
  public ClusterLabel(int, float, int, boolean);
  public java.util.Map<java.lang.Integer, java.lang.Float> getLabel();
  public float getWeight();
  public java.util.TreeSet<java.lang.Integer> labelSet();
  public java.util.Map<java.lang.Integer, java.lang.Float> labelMap();
  public void add(int);
  public boolean sameAs(ClusterLabel);
  public void neighbour(ClusterLabel, float);
  public void noMore();
  public void normalize();
  public java.lang.String toString();
}
Compiled from "COPRA.java"
class LabelPair {
  public float value;
  public int pos;
  public LabelPair(float, int);
  public java.lang.String toString();
}
Compiled from "ModOverlap.java"
public class ModOverlap {
  static int maxVertex;
  public ModOverlap();
  public static void main(java.lang.String[]);
  public static double modOverlap(java.lang.String, java.lang.String);
  public static java.util.HashMap<java.lang.String, java.util.HashSet<java.lang.String>> readGraphEdges(java.lang.String);
  public static void normalize(java.util.HashMap<java.lang.Integer, java.util.HashMap<java.lang.String, java.lang.Double>>, java.util.HashMap<java.lang.String, java.lang.Double>);
  public static double modOverlap(java.util.HashMap<java.lang.Integer, java.util.HashMap<java.lang.String, java.lang.Double>>, java.util.HashMap<java.lang.String, java.util.HashSet<java.lang.String>>, double, double);
  static {};
}
Compiled from "Project.java"
public class Project {
  public Project();
  public static void main(java.lang.String[]);
}
Compiled from "COPRA.java"
class SetPair implements java.lang.Comparable<SetPair> {
  public int id;
  public java.util.TreeSet<java.lang.Integer> set;
  public SetPair(int, java.util.TreeSet<java.lang.Integer>);
  public int compareTo(SetPair);
  public int compareTo(java.lang.Object);
}
Compiled from "COPRA.java"
class VecPair {
  public java.util.Vector<java.lang.String> name;
  public java.util.Vector<java.lang.Double> value;
  public VecPair(java.util.Vector<java.lang.String>, java.util.Vector<java.lang.Double>);
}
Compiled from "COPRA.java"
class Vert implements java.lang.Comparable<Vert> {
  public int id;
  public int degree;
  public Vert(int, int);
  public int compareTo(Vert);
  public int compareTo(java.lang.Object);
}
