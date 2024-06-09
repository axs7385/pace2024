#include <bits/stdc++.h>
#include <signal.h>
#include <unistd.h>

#define sig
using namespace std;
#define ll long long
#ifdef sig
volatile sig_atomic_t tle = 0;
void term(int signum)
{
    tle = 1;
}
#else
int tle = 0;
#endif
const int mod = 998244353, base = 19260817; // hash å¸¸æ°

set<ll> tabu;                    // ç¦å¿è¡¨
const int N = 300005, M = 21000; // æ»åº
int n, m, K;
vector<int> e[N];
// int c[M][M],p[M+1][M+1];
int ord[N], _ord[N];
int ords[10][M + 5], used[M + 5], child[M + 5];
ll ans[10], childans;
vector<ll> c[M], p[M + 1];
double avg[N];

// ååºæå
bool cmp(int x, int y)
{
    return avg[x] < avg[y];
}

// è¯»å¥æ°æ®
void read(int &x)
{
    x = 0;
    char c = getchar();
    while (c == 'c')
    {
        while (c != '\n')
            c = getchar();
        c = getchar();
    }
    while (c > '9' || c < '0')
        c = getchar();
    while (c >= '0' && c <= '9')
    {
        x = x * 10 + c - '0';
        c = getchar();
    }
}

ll calc(int i, int j)
{
    int l = 0;
    int lj = e[j].size();
    ll sum = 0;
    for (int x : e[i])
    {
        while (l < lj && e[j][l] < x)
            ++l;
        sum += l;
    }
    return sum;
}

ll get_ans(int ord[])
{
    ll sum = 0;
    for (int i = 1; i <= m; ++i)
        for (int j = i + 1; j <= m; ++j)
            sum += c[ord[i]][ord[j]];
    return sum;
}

ll dp[N];
int fr[N], opt[N];

//å¯¹ordè¿ä¸ªè§£è¿è¡ä¸éæç´¢ï¼åºäºPengçdpæ¹æ¡
ll work(int ord[],int bg=0)
{
    for (int i = 1; i <= m; ++i)
    {
        for (int j = 1; j <= m; ++j)
            p[i][j] = p[i][j - 1] + c[ord[i]][ord[j]];
    }
    dp[0] = 0;
    for (int i=1;i<bg;++i)
    {
        dp[i]=dp[i-1];
        fr[i]=1;
        opt[i]=0;
    }
    for (int i = bg; i <= m; ++i)
    {
        dp[i] = dp[i - 1];
        fr[i] = 1;
        opt[i] = 0;
        // cerr<<"."<<i;
        // if (i<=bg)
        //     continue;
        for (int j = 1; j <= i - 1; ++j)
        {
            //å°ç¬¬jä¸ªç¹æå¥å°ç¬¬iä¸ªç¹çåä¸ä¸ªä½ç½®
            if (dp[j - 1] + (p[j][i] - p[j][j]) > dp[i])
            {
                dp[i] = dp[j - 1] + (p[j][i] - p[j][j]);
                fr[i] = i - j + 1;
                opt[i] = 1;
            }
            //å°ç¬¬iä¸ªç¹æå¥å°ç¬¬jä¸ªç¹çåä¸ä¸ªä½ç½®
            if (dp[j - 1] - (p[i][i] - p[i][j - 1]) > dp[i])
            {
                dp[i] = dp[j - 1] - (p[i][i] - p[i][j - 1]);
                fr[i] = i - j + 1;
                opt[i] = 2;
            }
            ////å°ç¬¬iä¸ªç¹åç¬¬jä¸ªç¹äº¤æ¢
            if (dp[j - 1] + (p[j][i] - p[j][j]) - (p[i][i] - p[i][j]) > dp[i])
            {
                dp[i] = dp[j - 1] + (p[j][i] - p[j][j]) - (p[i][i] - p[i][j]);
                fr[i] = i - j + 1;
                opt[i] = 3;
            }

            if (j+3<i)
            {
                 //å°jåj+1æå¥å°iå
                if (dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1])>dp[i])
                {
                    dp[i]=dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1]);
                    fr[i]=i-j+1;
                    opt[i]=5;
                }
                 //å°i-1åiæå¥å°jå
                if (dp[j-1]-(p[i][i-2]-p[i][j-1])-(p[i-1][i-2]-p[i-1][j-1])>dp[i])
                {
                    dp[i]=dp[j-1]-(p[i][i-2]-p[i][j-1])-(p[i-1][i-2]-p[i-1][j-1]);
                    fr[i]=i-j+1;
                    opt[i]=6;
                }
            }
            if (j+4<i)
            {
                //å°j,j+1åi-1,iäº¤æ¢
                if (dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1])-(p[i][i-2]-p[i][j+1])-(p[i-1][i-2]-p[i-1][j+1])>dp[i])
                {
                    dp[i]=dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1])-(p[i][i-2]-p[i][j+1])-(p[i-1][i-2]-p[i-1][j+1]);
                    fr[i]=i-j+1;
                    opt[i]=7;
                }
                //å°j,j+1åiäº¤æ¢
                if (dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1])-(p[i][i]-p[i][j+1])>dp[i])
                {
                    dp[i]=dp[j-1]+(p[j][i]-p[j][j+1])+(p[j+1][i]-p[j+1][j+1])-(p[i][i]-p[i][j+1]);
                    fr[i]=i-j+1;
                    opt[i]=8;
                }
                //å°jåi-1,iäº¤æ¢
                if (dp[j-1]+(p[j][i]-p[j][j])-(p[i][i-2]-p[i][j])-(p[i-1][i-1]-p[i-1][j])>dp[i])
                {
                    dp[i]=dp[j-1]+(p[j][i]-p[j][j])-(p[i][i-2]-p[i][j])-(p[i-1][i-1]-p[i-1][j]);
                    fr[i]=i-j+1;
                    opt[i]=9;
                }
            }
            
        }
    }
    
    if (dp[m] <= 0)
        return 0;
    int cnt = m;
    for (int i = m; i; i -= fr[i])
    {
        // cerr<<"?"<<i;
        if (fr[i] == 1)
        {
            _ord[cnt--] = ord[i];
            continue;
        }
        else
        {
            // cerr<<'('<<i-fr[i]+1<<' '<<i<<' '<<opt[i]<<' '<<dp[i-fr[i]]<<' '<<dp[i]<<')'<<endl;
            if (opt[i]<=4)
            {
                if (opt[i] & 1)
                    _ord[cnt--] = ord[i - fr[i] + 1];
                if ((opt[i] & 2) == 0)
                    _ord[cnt--] = ord[i];
                for (int j = i - 1; j > i - fr[i] + 1; --j)
                    _ord[cnt--] = ord[j];
                if ((opt[i] & 1) == 0)
                    _ord[cnt--] = ord[i - fr[i] + 1];
                if (opt[i] & 2)
                    _ord[cnt--] = ord[i];
            }
            if (opt[i]==5)
            {
                _ord[cnt--]=ord[i-fr[i]+2];
                _ord[cnt--]=ord[i-fr[i]+1];
                for (int j=i;j>i-fr[i]+2;--j)
                    _ord[cnt--]=ord[j];
            }
            if (opt[i]==6)
            {
                for (int j=i-2;j>i-fr[i];--j)
                    _ord[cnt--]=ord[j];
                _ord[cnt--] = ord[i];
                _ord[cnt--] = ord[i-1];
            }
            if (opt[i]==7)
            {
                _ord[cnt--]=ord[i-fr[i]+2];
                _ord[cnt--]=ord[i-fr[i]+1];
                for (int j=i-2;j>i-fr[i]+2;--j)
                    _ord[cnt--]=ord[j];
                _ord[cnt--] = ord[i];
                _ord[cnt--] = ord[i-1];
            }
            if (opt[i]==8)
            {
                _ord[cnt--]=ord[i-fr[i]+2];
                _ord[cnt--]=ord[i-fr[i]+1];
                for (int j=i-1;j>i-fr[i]+2;--j)
                    _ord[cnt--]=ord[j];
                _ord[cnt--] = ord[i];
                // _ord[cnt--] = ord[i-1];
            }
            if (opt[i]==9)
            {
                // _ord[cnt--]=ord[i-fr[i]+2];
                _ord[cnt--]=ord[i-fr[i]+1];
                for (int j=i-2;j>i-fr[i]+1;--j)
                    _ord[cnt--]=ord[j];
                _ord[cnt--] = ord[i];
                _ord[cnt--] = ord[i-1];
            }
        }
    }
    // cerr<<"!"<<cnt<<endl;
    for (int i = 1; i <= m; ++i)
        ord[i] = _ord[i];
    // cerr<<"!"<<dp[m]<<endl;
    return dp[m];
}
ll workbig(int ord[])
{
    dp[0] = 0;
    for (int i = 1; i <= m; ++i)
        dp[i] = -1;
    for (int i = 1; i <= m&&!tle; ++i)
    {
        // if (tle)
        //     return 0;
        // cout<<i<<endl;
        if (dp[i] < dp[i - 1])
        {
            opt[i] = 0;
            fr[i] = 1;
            dp[i] = dp[i - 1];
        }
        ll sum = 0;
        for (int j = i - 1; j&&!tle; --j)
        {
            sum -= calc(ord[i], ord[j]) - calc(ord[j], ord[i]);
            if (dp[j - 1] + sum > dp[i])
            {
                dp[i] = dp[j - 1] + sum;
                fr[i] = i - j + 1;
                opt[i] = 2;
            }
        }
        sum = 0;
        for (int j = i + 1; j <= m&&!tle; ++j)
        {
            sum += calc(ord[i], ord[j]) - calc(ord[j], ord[i]);
            if (dp[i - 1] + sum > dp[j])
            {
                dp[j] = dp[i - 1] + sum;
                fr[j] = j - i + 1;
                opt[j] = 1;
            }
        }
    }
    if (tle)
        return 0;
    if (dp[m] <= 0)
        return 0;
    int cnt = m;
    for (int i = m; i; i -= fr[i])
    {
        if (fr[i] == 1)
        {
            _ord[cnt--] = ord[i];
            continue;
        }
        else
        {
            if (opt[i] & 1)
                _ord[cnt--] = ord[i - fr[i] + 1];
            if ((opt[i] & 2) == 0)
                _ord[cnt--] = ord[i];
            for (int j = i - 1; j > i - fr[i] + 1; --j)
                _ord[cnt--] = ord[j];
            if ((opt[i] & 1) == 0)
                _ord[cnt--] = ord[i - fr[i] + 1];
            if (opt[i] & 2)
                _ord[cnt--] = ord[i];
        }
    }
    for (int i = 1; i <= m; ++i)
        ord[i] = _ord[i];
    // cout<<dp[m]<<':';for (int i=1;i<=m;++i)cout<<ord[i]<<' ';cout<<endl;
    return dp[m];
}
mt19937 rnd(0);
const int nn = 18;
ll f[1 << nn];
int _fr[1 << nn];
const ll inf = 1e15;
int changeord(int ord[])
{

    int len = min(nn, m);
    int x;
    if (len == m)
        x = 0;
    else
        x = rnd() % (m - len);
    x += 1;
    ll sum = 0;
    for (int i = 0; i < len; ++i)
        for (int j = i + 1; j < len; ++j)
            sum += c[ord[i + x]][ord[j + x]];
    f[0] = 0;
    for (int i = 1; i < (1 << len); ++i)
    {
        f[i] = inf;
        for (int j = 0; j < len; ++j)
            if ((i >> j) & 1)
            {
                ll det = 0;
                for (int k = 0; k < len; ++k)
                    if ((i >> k) & 1)
                        det += c[ord[k + x]][ord[j + x]];
                if (det + f[i ^ (1 << j)] < f[i])
                {
                    _fr[i] = j;
                    f[i] = f[i ^ (1 << j)] + det;
                }
            }
    }
    if (f[(1 << len) - 1] == sum)
        return 0;
    vector<int> nord;
    for (int i = (1 << len) - 1; i;)
    {
        nord.push_back(ord[_fr[i] + x]);
        i -= (1 << _fr[i]);
    }
    assert(nord.size() == len);
    for (int i = 0; i < len; ++i)
        ord[i + x] = nord[len - i - 1];

    // set<int> ss;
    // for (int i=0;i<n)
    return x;
}
int swapsum=0;
void cross(int parent1[], int parent2[], int child[])
{
    int x = rnd() % m + 1;
    int y = rnd() % m + 1;
    if (x > y)
        swap(x, y);
    for (int i = 0; i < m; ++i)
        used[i] = 0;
    for (int i = x; i <= y; ++i)
    {
        child[i] = parent1[i];
        used[child[i]] = 1;
    }
    int pos = 0;
    for (int i = 1; i <= m; ++i)
    {
        if (!used[parent2[i]])
        {
            ++pos;
            if (pos == x)
                pos = y + 1;
            child[pos] = parent2[i];
        }
    }
    // for (int i=1;i<=swapsum;++i)
    // {
    //     int x=rnd()%m+1;int y=rnd()%m+1;
    //     if (x!=y)
    //         swap(child[x],child[y]);
    // }
    // swapsum+=5;
    // if (swapsum>n/5)
    //     swapsum=0;
}
ll get_hash(int ord[])
{
    ll sum = 0;
    for (int i = 1; i <= m; ++i)
        sum = ((ll)sum * base % mod + ord[i]) % mod;
    return sum;
}

int main(int argc, char *argv[])
{
    // freopen("1.gr","r",stdin);
#ifdef sig
    struct sigaction action; // ä¿¡å·
    memset(&action, 0, sizeof(struct sigaction));
    action.sa_handler = term;
    sigaction(SIGTERM, &action, NULL);
#endif
    clock_t clockBegin = clock();
    read(n);                          // Aç¹æ°
    read(m);                          // Bç¹æ°
    read(K);                          // è¾¹æ°
    for (int i = 0, u, v; i < K; ++i) // æ¯ä¸ªç¹ä»0å¼å§ç¼å·
    {
        read(u);
        --u;
        read(v);
        --v;
        v -= n;
        e[v].push_back(u); // eä¿å­äºææçè¾¹ï¼e[v]è¡¨ç¤ºvçé»å± v in B
    }
    for (int i = 0; i < m; ++i) // éåb
        if (!e[i].empty())
        {
            sort(e[i].begin(), e[i].end()); // içé»å±ååºæåº
            for (int x : e[i])
                avg[i] += x;
            avg[i] /= (double)e[i].size(); // avg[i]æ¯içé»å±çå¹³åæ°
        }
        else
        {
            avg[i] = 0;
            // cerr<<"!"<<i<<endl;
        }
    for (int i = 1; i <= m; ++i)
        ord[i] = i - 1;
    sort(ord + 1, ord + 1 + m, cmp); // Bçé¡ºåºï¼åå§è§£ï¼ cmpæ¯è¾çæ¯avg
    // cerr<<m<<' '<<M<<endl;
    /**
    å°B>21000æ¶ï¼è®¤ä¸ºåå­å­ä¸ä¸ï¼N^2ç©ºé´
    */
    if (m < M)
    // if (0)
    {
        for (int i = 0; i < m; ++i)
        {
            c[i].resize(m);
            p[i].resize(m + 1);
        }
        p[m].resize(m + 1);
        for (int i = 0; i < m; ++i)
            for (int j = i + 1; j < m; ++j)
            {
                ll pi = calc(i, j), pj = calc(j, i);
                c[i][j] = pi - pj; // C[i][j] ç¹iå¨jä¹åçï¼N[i]åN[j]äº§ççäº¤åæ°
                c[j][i] = pj - pi; //
            }
        // cerr<<"finish"<<' '<<clock()-clockBegin<<endl;
        int mxt = 290;               // 29*60; //æ¶é´cutoff
        // cerr<<60*30*CLOCKS_PER_SEC<<endl;
        int mx = 0;                  // ç§ç¾¤å¤§å°ï¼ä¸æ¯å¸¸æ°ï¼ä¼ååï¼æå¤§10
        for (int i=0;i==0||(i<10&&clock()-clockBegin<=mxt/5*CLOCKS_PER_SEC);++i) // æå¤§çç§ç¾¤ä¸º10ï¼æé 10æ¬¡åå§è§£
        {
            //  cerr << i << endl;
            for (int j = 1; j <= m; ++j)
                ords[i][j] = ord[j]; // orders[i]è¡¨ç¤ºç¬¬iä¸ªä¸ªä½
            if (i >= 1) //é¤äºorders[0]ï¼å¶ä»çè§£é½éæºå
            {
                // shuffle(ords[i] + 1, ords[i] + 1 + m, rnd);
                for (int j=1;j<=m;++j)
                    swap(ords[i][j],ords[i][j+rnd()%min(m-j+1,i*10)]);
            }

            ll x = get_hash(ords[i]); 
            int cnt = 0;
            /**
             * tabu.count(x) != 0 ->xåºç°è¿            
            */
            while (tabu.count(x) != 0 && cnt <= 10)// && !tle && clock() - clockBegin <= mxt / 2 * CLOCKS_PER_SEC)
            {
                shuffle(ords[i] + 1, ords[i] + 1 + m, rnd);
                x = get_hash(ords[i]);
                ++cnt;
            }
            // cerr<<"hashb"<<x<<endl;
            //ç»ææé åå§ç§ç¾¤
            if (tabu.count(x) != 0 || tle)// || clock() - clockBegin > mxt / 2 * CLOCKS_PER_SEC)
                break;
            // tabu.insert(x);
            // ll tmp;
            // while ((tmp=work(ords[i]))&& !tle)// && clock() - clockBegin <= mxt / 2 * CLOCKS_PER_SEC)
            //     cerr<<tmp<<' '<<(clock()-clockBegin)/(double)CLOCKS_PER_SEC<<';';

            while (work(ords[i]) && !tle)
                ;
            x = get_hash(ords[i]);
            while (tabu.count(x) != 0 && cnt <= 10 && !tle)// && !tle & clock() - clockBegin <= mxt / 2 * CLOCKS_PER_SEC)
            {
                shuffle(ords[i] + 1, ords[i] + 1 + m, rnd);
                while (work(ords[i]) && !tle)// && clock() - clockBegin <= mxt / 2 * CLOCKS_PER_SEC)
                    ;
                x = get_hash(ords[i]);
                ++cnt;
            }
            if (tle)
                break;
            if (tabu.count(x))
                break;
            
            tabu.insert(x);
            if (!tle)
            {
                mx = i + 1;
                ans[i] = get_ans(ords[i]);
            }
            // cerr<<"hash:"<<x<<' '<<ans[i]<<endl;
            // cerr<<(clock()-clockBegin)/(double)CLOCKS_PER_SEC<<endl;
        }
        // cerr<<"!"<<clock()-clockBegin<<endl;
        if (mx == 0)
        {
            for (int i = 1; i <= m; ++i)
                printf("%d\n", ord[i] + n + 1);
            return 0;
        }
        int pc=0;
        while (clock() - clockBegin <= mxt / 2 * CLOCKS_PER_SEC && !tle && mx > 1)
        {
            // cerr<<(clock()-clockBegin)/(double)CLOCKS_PER_SEC<<endl;
            // cerr<<'.';
            int x = rnd() % mx, y = rnd() % mx;
            while (x == y)
                y = rnd() % mx;
            cross(ords[x], ords[y], child);
            int z = get_hash(child);
            if (tabu.count(z))
                continue;
            while (work(child) && !tle)// && clock() - clockBegin <= mxt / 2 * CLOCKS_PER_SEC)
                ;
            ll tmphash = get_hash(child);
            if (tabu.count(tmphash))
                continue;
            int mnp = 0;
            tabu.insert(tmphash);
            for (int i = 1; i < mx; ++i)
                if (ans[i] > ans[mnp])// || (ans[i] == ans[mnp] && rnd() % 2 == 0))
                    mnp = i;
            ll tmpans = get_ans(child);
            // if (++pc<=10)cerr<<tmphash<<' '<<tmpans<<endl;
            if (tmpans <= ans[mnp])
            {
                ans[mnp] = tmpans;
                for (int i = 1; i <= m; ++i)
                    ords[mnp][i] = child[i];
                // cout<<tmpans<<' '<<tmphash<<' ';for (int i=1;i<=m;++i)cout<<child[i]<<' ';cout<<endl;
            }
        }
        int mxp = 0;
        for (int i = 1; i < mx; ++i)
            if (ans[i] < ans[mxp])
                mxp = i;
        for (int i = 1; i <= m; ++i)
            ord[i] = ords[mxp][i];
        // cerr<<"!"<<get_ans(ord)<<endl;
        // while (clock() - clockBegin <= mxt * CLOCKS_PER_SEC && !tle)
        // {

        //     if (!work(ord))
        //     {
        //         while (clock() - clockBegin <= mxt * CLOCKS_PER_SEC && !tle)
        //         {
        //             if (changeord(ord))
        //                 break;
        //         }
        //     }
        //     // break;
        // }
        int bg=0;
        while (clock() - clockBegin <= mxt * CLOCKS_PER_SEC && !tle)
        {
            // cerr<<clock()-clockBegin<<endl;
            if (!work(ord,bg))
            {
                while (clock() - clockBegin <= mxt * CLOCKS_PER_SEC && !tle)
                {
                    int pos=changeord(ord);
                    if (pos!=0)
                    {
                        bg=pos;
                        break;
                    }
                }
            }
            else
                bg=0;
            // break;
        }
    }
    else
    {
        while (clock() - clockBegin <= 5 * CLOCKS_PER_SEC && !tle)
        {
            if (!workbig(ord))
                break;
        }
    }
    for (int i = 1; i <= m; ++i)
        printf("%d\n", ord[i] + n + 1);
    // printf("%lld\n",get_ans(ord));
    // cerr<<"!"<<get_ans(ord)<<endl;
    return 0;
}