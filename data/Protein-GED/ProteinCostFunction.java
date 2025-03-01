/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package util;

/**
 *
 * @author mferrer
 */

import util.GraphComponent;
import algorithms.Constants;
import algorithms.MunkresRec;

public class ProteinCostFunction implements ICostFunction {

	private double[][] stringMatrix;

	/**
	 * the constant costs
	 */
	private double nodeCosts;

	private double edgeCosts;
	
	private double alpha;

	private MunkresRec munkresRec;

	public ProteinCostFunction(double nodeCosts, double edgeCosts, double alpha) {
		this.nodeCosts = nodeCosts;
		this.edgeCosts = edgeCosts;
		this.alpha = alpha;
		this.munkresRec = new MunkresRec();
	}

	/**
	 * @return costs of a distortion between
	 * @param start
	 *            and
	 * @param end
	 */
	public double getCosts(GraphComponent start, GraphComponent end) {
		
		if( Constants.nodecostmatrix!=null && Constants.edgecostmatrix!=null){
			return precomputedcosts(start,end);
		}
		
		/**
		 * node handling
		 */
		if (start.isNode() || end.isNode()) {
			String chemSym1;
			String chemSym2;
			if (!start.getComponentId().equals(Constants.EPS_ID)) { // start is not empty
				chemSym1 = (String) start.getTable().get("type");
			} else { // insertion
				return this.alpha * 1 * this.nodeCosts;
			}

			if (!end.getComponentId().equals(Constants.EPS_ID)) {// end is not empty
				chemSym2 = (String) end.getTable().get("type");
			} else { // deletion
				return this.alpha * 1* this.nodeCosts;
			}
			
			if(chemSym1.equals(chemSym2)==true){
				String seq1 = (String) start.getTable().get("sequence");
				String seq2 = (String) end.getTable().get("sequence");
				return this.alpha * this.getStringDistance(seq1, seq2);
			}else{
				return this.alpha * 1 * this.nodeCosts;
			}
		
		}
		/**
		 * edge handling
		 */
		else {
			int startFreq;
			int endFreq;
			if (start.getComponentId().equals(Constants.EPS_ID)) {
				endFreq = Integer.parseInt((String) end.getTable().get("frequency"));
				return (1-this.alpha)*endFreq* this.edgeCosts;
			}
			if (end.getComponentId().equals(Constants.EPS_ID)) {
				startFreq = Integer.parseInt((String) start.getTable().get("frequency"));
				return (1-this.alpha)*startFreq * this.edgeCosts;
			}
			startFreq = Integer.parseInt((String) start.getTable().get("frequency"));
			endFreq = Integer.parseInt((String) end.getTable().get("frequency"));
			int n = startFreq + endFreq;
			double[][] matrix = new double[n][n];
			for (int i = 0; i < startFreq; i++){
				String startType = (String) start.getTable().get("type"+i);
				for (int j = 0; j < endFreq; j++){
					String endType = (String) end.getTable().get("type"+j);
					if (startType.equals(endType)){
						matrix[i][j] = 0;
					} else {
						matrix[i][j] = 2*this.edgeCosts;
					}
				}
			}
			for (int i = startFreq; i < startFreq+endFreq; i++){
				for (int j = 0; j < endFreq; j++){
					if (i-startFreq == j){
						matrix[i][j] = this.edgeCosts;
					} else {
						matrix[i][j] = Double.POSITIVE_INFINITY;
					}
				}
			}
			for (int i = 0; i < startFreq; i++){
				for (int j = endFreq; j < startFreq + endFreq; j++){
					if (j-endFreq == i){
						matrix[i][j] = this.edgeCosts;
					} else {
						matrix[i][j] = Double.POSITIVE_INFINITY;
					}
				}
			}
			return (1-this.alpha)*this.munkresRec.getCosts(matrix);

			
		}
	}

	private double getStringDistance(String s1, String s2) {
		int n = s1.length();
		int m = s2.length();
		if (m > n) {
			String s = s1;
			s1 = s2;
			s2 = s;
			n = s1.length();
			m = s2.length();
		}
		s2 += s2;
		m *= 2.;
		this.stringMatrix = new double[n + 1][m + 1];
		this.stringMatrix[0][0] = 0;
		for (int i = 1; i <= n; i++) {
			this.stringMatrix[i][0] = this.stringMatrix[i - 1][0]
					+ this.nodeCosts;
		}
		for (int j = 1; j <= m; j++) {
			this.stringMatrix[0][j] = this.stringMatrix[0][j - 1];
		}

		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= m; j++) {
				double subst = 0.;
				if (s1.charAt(i - 1) == s2.charAt(j - 1)) {
					subst = 0.;
				} else {
					subst = this.nodeCosts;
				}
				double m1 = this.stringMatrix[i - 1][j - 1] + subst;
				double m2 = this.stringMatrix[i - 1][j] + this.nodeCosts;
				double m3 = this.stringMatrix[i][j - 1] + this.nodeCosts;
				this.stringMatrix[i][j] = Math.min(m1, Math.min(m2, m3));
			}
		}
		double min = Double.POSITIVE_INFINITY;
		for (int j = 0; j <= m; j++) {
			double current = this.stringMatrix[n][j];
			if (current < min) {
				min = current;
			}
		}
		return min;
	}

	/**
	 * @return the cost of an edge operation
	 */
	public double getEdgeCosts() {
		return (1-this.alpha) * edgeCosts;
	}

	/**
	 * @return the cost of a node operation
	 */
	public double getNodeCosts() {
		return this.alpha * nodeCosts;
	}
	
	private double precomputedcosts(GraphComponent start, GraphComponent end) {
		// TODO Auto-generated method stub
		
		double[][] matrix=null;
		
		if (start.isNode() || end.isNode()) {
			matrix=Constants.nodecostmatrix;
		}else{
			matrix=Constants.edgecostmatrix;
		}
		
		int n1 = matrix.length;
		int n2 = matrix[0].length;
		int insertindexg1=n2-2;
		int insertindexg2=n1-2;
		int delindexg1=n2-1;
		int delindexg2=n1-1;


		if (start.getComponentId().equals(Constants.EPS_ID)) {
			if(end.belongtosourcegraph){
				return matrix[end.id][insertindexg1];
			}else{
				return matrix[insertindexg2][end.id];
			}
		}

		if (end.getComponentId().equals(Constants.EPS_ID)) {
			if(start.belongtosourcegraph){
				return matrix[start.id][delindexg1];
			}else{
				return matrix[delindexg2][start.id];
			}
		}
		
		if(start.belongtosourcegraph){
			return matrix[start.id][end.id];
		}else{
			return matrix[end.id][start.id];
		}




	}

}

