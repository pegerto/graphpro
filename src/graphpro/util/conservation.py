class ConservationScoreClient:
    def __init__(self, pdb_id: str, chain_id: str, server: str = 'graphpro.pegerto.com'):
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.server = server

    def compute_conservation_score(self) -> list[tuple[int,float]]:
        import requests
        
        resouce = f'https://{self.server}/conservation/{self.pdb_id}/{self.chain_id}'
        resp = requests.post(resouce, verify=False)
        if resp.status_code != 200:
            raise Exception(f"Error calculating conservation: {resp.status_code}")
        
        results = resp.json()['chains'][self.chain_id]
        return [(int(resid), float(score['shannon'])) for resid, score in results.items()]