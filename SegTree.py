class Seg_min():
    """
    0-indexedの配列[a0, a1, a2, ..., a(x-1)]に対して以下のクエリをそれぞれO(logx)で行う：
        1. i番目の要素にxを代入
        2. 区間[i, j]内の最小値を返す
    """
    def __init__(self,x):
        #####単位元######
        self.ide_ele_min = 10**10
        self.func = min

        self.n = len(x)

        #num_max:n以上の最小の2のべき乗
        self.num_max =2**(self.n-1).bit_length()
        self.x = [self.ide_ele_min]*2*self.num_max

        for i,num in enumerate(x, self.num_max):
            self.x[i] = num
        for i in range(self.num_max-1,0,-1):
            self.x[i] = self.func(self.x[i<<1],self.x[(i<<1) + 1])

    def update(self,i,x):
        """
        i番目の要素にxを代入
        """
        i += self.num_max
        self.x[i] = x
        while(i>0):
            i = i//2
            self.x[i] = self.func(self.x[i<<1],self.x[(i<<1) + 1])

    def query(self,i,j):
        """
        区間[i, j]内の最小値を返す
        """
        res = self.ide_ele_min
        if i>=j:
            return res
        i += self.num_max
        j += self.num_max -1
        while(i<=j):
            if(i==j):
                res = self.func(res,self.x[i])
                break
            if(i&1):
                res = self.func(res,self.x[i])
                i += 1
            if(not j&1):
                res = self.func(res,self.x[j])
                j -= 1
            i = i>>1
            j = j>>1
        return res
