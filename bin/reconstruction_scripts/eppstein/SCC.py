Path: main.gmane.org!not-for-mail
From: "Tim Peters" <tim_one@email.msn.com>
Newsgroups: gmane.comp.python.devel
Subject: RE: Garbage collecting closures
Date: Wed, 16 Apr 2003 00:51:27 -0400
Lines: 185
Sender: python-dev-admin@python.org
Approved: news@gmane.org
Message-ID: <LNBBLJKPBEHFEDALKOLCEEBCEHAB.tim_one@email.msn.com>
References: <200304151923.h3FJNmG29436@odiug.zope.com>
NNTP-Posting-Host: main.gmane.org
Mime-Version: 1.0
Content-Type: multipart/mixed;
	boundary="----=_NextPart_000_0006_01C303B2.4FA6CF60"
X-Trace: main.gmane.org 1050468747 27594 80.91.224.249 (16 Apr 2003 04:52:27 GMT)
X-Complaints-To: usenet@main.gmane.org
NNTP-Posting-Date: Wed, 16 Apr 2003 04:52:27 +0000 (UTC)
Original-X-From: python-dev-admin@python.org Wed Apr 16 06:52:25 2003
Return-path: <python-dev-admin@python.org>
Original-Received: from mail.python.org ([12.155.117.29])
	by main.gmane.org with esmtp (Exim 3.35 #1 (Debian))
	id 195euH-0007Ao-00
	for <python-python-dev@m.gmane.org>; Wed, 16 Apr 2003 06:52:25 +0200
Original-Received: from localhost.localdomain ([127.0.0.1] helo=mail.python.org)
	by mail.python.org with esmtp (Exim 4.05)
	id 195euw-0000Az-00; Wed, 16 Apr 2003 00:53:06 -0400
Original-Received: from bay0-smtp09.bay0.hotmail.com ([65.54.241.116] helo=BAY0-SMTP09.adinternal.hotmail.com)
	by mail.python.org with esmtp (Exim 4.05)
	id 195ett-00009S-00
	for python-dev@python.org; Wed, 16 Apr 2003 00:52:01 -0400
X-Originating-IP: [204.30.139.38]
X-Originating-Email: [tim_one@msn.com]
Original-Received: from cj569191b ([204.30.139.38]) by BAY0-SMTP09.adinternal.hotmail.com with Microsoft SMTPSVC(5.0.2195.5600);
	 Tue, 15 Apr 2003 21:51:28 -0700
Original-To: <python-dev@python.org>
X-Priority: 3 (Normal)
X-MSMail-Priority: Normal
X-Mailer: Microsoft Outlook IMO, Build 9.0.2416 (9.0.2911.0)
In-reply-to: <200304151923.h3FJNmG29436@odiug.zope.com>
X-MIMEOLE: Produced By Microsoft MimeOLE V6.00.2800.1106
Importance: Normal
X-OriginalArrivalTime: 16 Apr 2003 04:51:29.0121 (UTC) FILETIME=[D7DF1510:01C303D3]
Errors-To: python-dev-admin@python.org
X-BeenThere: python-dev@python.org
X-Mailman-Version: 2.0.13 (101270)
Precedence: bulk
List-Help: <mailto:python-dev-request@python.org?subject=help>
List-Post: <mailto:python-dev@python.org>
List-Subscribe: <http://mail.python.org/mailman/listinfo/python-dev>,
	<mailto:python-dev-request@python.org?subject=subscribe>
List-Id: Python core developers <python-dev.python.org>
List-Unsubscribe: <http://mail.python.org/mailman/listinfo/python-dev>,
	<mailto:python-dev-request@python.org?subject=unsubscribe>
List-Archive: <http://mail.python.org/pipermail/python-dev/>
Xref: main.gmane.org gmane.comp.python.devel:11663
X-Report-Spam: http://spam.gmane.org/gmane.comp.python.devel:11663

[Guido]
> I'm glazing over the details now, but there seems to be a kernel of
> useful cleanup in here somehow; I hope that someone will be able to
> contribute a prototype of such code at least!

I'll attach a head start, a general implementation of Tarjan's SCC algorithm
that produces a list of SCCs already in a topsort order.  I haven't tested
this enough, and Tarjan's algorithm is subtle -- user beware.

The trygc() function at the end is an example application that appears to
work, busting all the objects gc knows about into SCCs and displaying them.
This requires Python CVS (for the new gc.get_referents function).  Note that
you'll get a very large SCC at the start.  This isn't an error!  Each module
that imports sys ends up in this SCC, due to that the module has the module
sys in its module dict, and sys has the module in its sys.modules dict.
From there, modules have their top-level functions in their dict, while the
top level functions point back to the module dict via func_globals.  Etc.
Everything in this giant blob is reachable from everything else.

For the gc application, it would probably be better (run faster and consume
less memory) if dfs() simply ignored objects with no successors.
Correctness shouldn't be harmed if def started with

    succs = successors(v)
    if not succs:
        return

except that objects with no successors would no longer be considered
singleton SCCs, and the recursive call to dfs() would need to be fiddled to
skip trying to update id2lowest[v_id] then (so dfs should be changed to
return a bool saying whether it took the early return).  This would save the
current work of trying to chase pointless things like ints and strings.
Still, it's pretty zippy as-is!
---------------------------------------------------------------------
# This implements Tarjan's linear-time algorithm for finding the maximal
# strongly connected components.  It takes time proportional to the sum
# of the number of nodes and arcs.
#
# Two functions must be passed to the constructor:
#     node2id     graph node -> a unique integer
#     successors  graph node -> sequence of immediate successor graph nodes
#
# Call method getsccs() with an iterable producing the root nodes of the graph.
# The result is a list of SCCs, each of which is a list of graph nodes.
# This is a partitioning of all graph nodes reachable from the roots,
# where each SCC is a maximal subset such that each node in an SCC is
# reachable from all other nodes in the SCC.  Note that the derived graph
# where each SCC is a single "supernode" is necessarily acyclic (else if
# SCC1 and SCC2 were in a cycle, each node in SCC1 would be reachable from
# each node in SCC1 and SCC2, contradicting that SCC1 is a maximal subset).
# The list of SCCs returned by getsccs() is in a topological sort order wrt
# this derived DAG.

class SCC(object):
    def   init  (self, node2id, successors):
        self.node2id = node2id
        self.successors = successors

    def getsccs(self, roots):
        import sys

        node2id, successors = self.node2id, self.successors
        get dfsnum = iter(xrange(sys.maxint)).next
        id2dfsnum = {}
        id2lowest = {}
        stack = []
        id2stacki = {}
        sccs = []

        def dfs(v, v id):
            id2dfsnum[v id] = id2lowest[v id] = v dfsnum = get dfsnum()
            id2stacki[v id] = len(stack)
            stack.append((v, v id))
            for w in successors(v):
                w id = node2id(w)
                if w id not in id2dfsnum:   # first time we saw w
                    dfs(w, w id)
                    id2lowest[v id] = min(id2lowest[v id], id2lowest[w id])
                else:
                    w dfsnum = id2dfsnum[w id]
                    if w dfsnum < v dfsnum and w id in id2stacki:
                        id2lowest[v id] = min(id2lowest[v id], w dfsnum)

            if id2lowest[v id] == v dfsnum:
                i = id2stacki[v id]
                scc = []
                for w, w id in stack[i:]:
                    del id2stacki[w id]
                    scc.append(w)
                del stack[i:]
                sccs.append(scc)

        for v in roots:
            v id = node2id(v)
            if v id not in id2dfsnum:
                dfs(v, v id)
        sccs.reverse()
        return sccs

 basic tests = """
>>> succs = {1: [2], 2: []}
>>> s = SCC(int, lambda i: succs[i])

The order in which the roots are listed doesn't matter:  we get the unique
topsort regardless.

>>> s.getsccs([1])
[[1], [2]]
>>> s.getsccs([1, 2])
[[1], [2]]
>>> s.getsccs([2, 1])
[[1], [2]]

But note that 1 isn't reachable from 2, so giving 2 as the only root won't
find 1.

>>> s.getsccs([2])
[[2]]

>>> succs = {1: [2],
...          2: [3, 5],
...          3: [2, 4],
...          4: [3],
...          5: [2]}
>>> s = SCC(int, lambda i: succs[i])
>>> s.getsccs([1])
[[1], [2, 3, 4, 5]]
>>> s.getsccs(range(1, 6))
[[1], [2, 3, 4, 5]]

Break the link from 4 back to 2.
>>> succs[4] = []
>>> s.getsccs([1])
[[1], [2, 3, 5], [4]]
"""

  test   = {'basic':  basic tests}

def  test():
    import doctest
    doctest.testmod()

if   name   == '  main  ':
     test()

def trygc():
    import gc
    gc.collect()
    s = SCC(id, gc.get referents)
    for scc in s.getsccs(gc.get objects()):
        if len(scc) == 1:
            continue
        print "SCC w/", len(scc), "objects"
        for x in scc:
            print "   ", hex(id(x)), type(x),
            if hasattr(x, "  name  "):
                print x.  name  ,
            print
