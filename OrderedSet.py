import sqlite3

"""　MultiSetも同様に実装できそうなのでいつか試したい。。。　"""

class OrderedSet():
    """
    順序付き集合のPython実装。次の操作をlogオーダーで実行する：
        - 集合 S に x を追加
        - S に属する x より大きい要素のうち最小の要素を取得
        - S に属する x より小さい要素のうち最大の要素を取得

    ナイーブ実装だとTLEするのでsqliteを使う裏技で回避してる。
    ここの実装だと整数か実数しか入れられないので、より複雑なことをするなら
        - 自前でカラムを増やしたテーブルを定義する
        - BITで座標圧縮してOrderedSetもどきを作る（要素数N<=10^6に限る）
        - クエリを逆読みしてUnionFindを使う（適用可能ケースは限られる）
    あたりを試すこと。

    参考：https://atcoder.jp/contests/abc217/submissions/25621717
    """
    def __init__(self, dtype:str="INT"):
        if dtype == "INT":
            self.dtype = int
        elif dtype == "REAL":
            self.dtype = float
        else:
            raise KeyError(f"dtype={dtype} is not supported")
        
        db = sqlite3.connect(':memory:', isolation_level=None)
        self.cur = db.cursor()
        # see https://www.sqlite.org/pragma.html for each pragma order
        # CRAETE INDEX is for faster search
        self.cur.executescript("""
            PRAGMA trusted_schema = OFF;
            PRAGMA journal_mode = OFF;
            PRAGMA synchronous = OFF;
            PRAGMA temp_store = memory;
            PRAGMA secure_delete = OFF;
            """)
        self.cur.execute(f"CREATE TABLE tbl(val {dtype})")   
        self.cur.execute("CREATE INDEX tbl_index ON tbl(val)")
    
    def size(self) -> int:
        """ Return the number of elements """
        self.cur.execute("select count(*) from tbl")
        return self.cur.fetchall()[0][0]
    
    def empty(self) -> bool:
        """ Check if the table is empty """
        return self.size() == 0
    
    def insert(self, x):
        """ Insert x into the table 'IF x IS NOT IN THERE' """
        x = self.dtype(x)  # convert numpy dtype to python dtype
        self.cur.execute("SELECT * FROM tbl WHERE val = ?", (x, ))
        if not self.cur.fetchall():
            self.cur.execute("INSERT INTO tbl VALUES(?)", (x, ))
    
    def erase(self, x):
        """ Erase x from the table """
        x = self.dtype(x)  # convert numpy dtype to python dtype
        self.cur.execute("DELETE FROM tbl where val = ?", (x, ))
    
    def lower_bound(self, x, equal:bool):
        """ Return the minimum value in the table, which is larger than (or equal to) x """
        x = self.dtype(x)  # convert numpy dtype to python dtype
        if equal:
            self.cur.execute("SELECT MIN(val) FROM tbl WHERE val >= ?", (x, ))
        else:
            self.cur.execute("SELECT MIN(val) FROM tbl WHERE val > ?", (x, ))
        return self.cur.fetchall()[0][0]
    
    def upper_bound(self, x, equal:bool):
        """ Return the maximum value in the table, which is smaller than (or equal to) x """
        x = self.dtype(x)  # convert numpy dtype to python dtype
        if equal:
            self.cur.execute("SELECT MAX(val) FROM tbl WHERE val <= ?", (x, ))
        else:
            self.cur.execute("SELECT MAX(val) FROM tbl WHERE val < ?", (x, ))
        return self.cur.fetchall()[0][0]
